#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2016 Junko Tsuji

import sys, os.path, fileinput
from optparse import OptionParser

# format strings
prestr = "%s_pre_%d"
matstr = "%s_%d"
bedline = "%s\t%d\t%d\t%s\n"

def trnaFastaInput(lines):
    name, seq = "", ""
    for x in lines:
        x = x.rstrip()
        if x.startswith(">"):
            if name != "" and seq != "":
                yield name, seq
            name = x[1:]
            seq = ""
        else:
            seq = x.upper()
    if name != "" and seq != "":
        yield name, seq

def loadPrecursor(fasta, loci, introns, pre_cnt, pre_seq):
    uniq = []
    leader, trailer = "", ""
    for name, seq in trnaFastaInput(fileinput.input(fasta)):
        trna, locus, intron, leader, trailer = name.split("|")
        if seq not in uniq:
            uniq.append(seq)
            cnt = pre_cnt.get(trna, 0) + 1
            pre_seq.setdefault(prestr % (trna,cnt), seq)
            pre_cnt[trna] = cnt
        sp = prestr % (trna, pre_cnt[trna])
        if intron != "no_intron":
            ib, ie = map(int, intron.split(":")[-1].split("-"))
            b = int(locus.split(":")[-1].split("-")[0])
            if introns.get(sp, []):
                raise Exception("can't handle more than 1 intron in a sequence")
            introns.update({sp: [ib - b, ie - b]})
        loci.update({locus: sp})
    leader = int(leader.split("=")[-1])
    trailer = int(trailer.split("=")[-1])
    return leader, trailer

def generateMature(leader, trailer, introns, pre_cnt, pre_seq,
                   mature_cnt, mature_seq, pre_to_mature):
    track = {}
    prebed, matbed = [], []
    trnas = sorted(pre_cnt.keys())
    for trna in trnas:
        his = "G" if trna.startswith("His-") else ""
        for i in xrange(pre_cnt[trna]):
            pid = prestr % (trna, i+1)
            seq = pre_seq[pid]
            L = len(seq)
            mature = seq[leader : -trailer]
            prebed.append(bedline % (pid, 0, leader, "5_leader"))
            prebed.append(bedline % (pid, L-trailer, L, "3_trailer"))
            ipos = -1
            if pid in introns:
                i_beg, i_end = introns[pid]
                mature = mature[:i_beg] + mature[i_end:]
                prebed.append(bedline % (pid, i_beg+leader, i_end+leader, "intron"))
                ipos = i_beg
            mature = his + mature + "CCA"
            if mature not in track:
                cnt = mature_cnt.get(trna, 0) + 1
                mid = matstr % (trna, cnt)
                mature_seq.setdefault(mid, mature)
                matbed.append(bedline % (mid, 0, len(mature), "mature_trna"))
                track.update({mature: mid})
                mature_cnt[trna] = cnt
                if ipos > -1:
                    i_beg = ipos - 1 + len(his)
                    i_end = ipos + 1 + len(his)
                    matbed.append(bedline % (mid, i_beg, i_end, "exon_junction"))
            mid = track[mature]
            pre_to_mature.update({pid: mid})
    return prebed, matbed

def reassignIds(bed, loci, pre_to_mature, fname):
    f = open(fname, "w")
    for x in open(bed):
        x = x.rstrip().split("\t")
        locus = "%s:%s-%s" % (x[0], x[1], x[2])
        pid = loci[locus]
        mid = pre_to_mature[pid]
        x[3] = "|".join([pid, mid])
        f.write("\t".join(x) + "\n")
    f.close()

def writeBed(bed, fname):
    f = open(fname, "w")
    for x in bed: f.write(x)
    f.close()

def writeFasta(cnt, seq, formatstr, fname):
    f = open(fname, "w")
    trnas = sorted(cnt.keys())
    for trna in trnas:
        for i in xrange(cnt[trna]):
            trna_id = formatstr % (trna, (i+1))
            f.write(">" + trna_id + "\n")
            f.write(seq[trna_id] + "\n")
    f.close()

def main(opts, args):
    # load data from tRNA precursor FASTA
    loci, introns, pre_cnt, pre_seq = {}, {}, {}, {}
    leader, trailer = loadPrecursor(args[1], loci, introns, pre_cnt, pre_seq)
    # process precursors to mature tRNAs
    mature_cnt, mature_seq, pre_to_mature = {}, {}, {}
    prebed, matbed = generateMature(leader, trailer, introns, pre_cnt, pre_seq, 
                                    mature_cnt, mature_seq, pre_to_mature)
    # reassign IDs for genomic tRNA annotation
    gbed = os.path.basename(args[0]).split(".")[0] + "_clean.bed"
    reassignIds(args[0], loci, pre_to_mature, gbed)

    # write BED files for tRNA reference
    writeBed(prebed, "precursor_trna.bed")
    writeBed(matbed, "mature_trna.bed")

    # write FASTA files for tRNA reference
    writeFasta(pre_cnt, pre_seq, prestr, "precursor_trna.fa")
    writeFasta(mature_cnt, mature_seq, matstr, "mature_trna.fa")

if __name__ == "__main__":
    usage = "%prog trnaGenomicBed trnaPrecursorFasta"
    description = "Prepare mature and precursor tRNA annotations and sequences"

    op = OptionParser(usage=usage, description=description)
    (opts, args) = op.parse_args()

    if len(args) != 2: op.error("input two arguments")

    try:
        main(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
