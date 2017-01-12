#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2016 Junko Tsuji

import sys, os.path, fileinput
from optparse import OptionParser

def inputReadStats(args):
    w, l = {}, set()
    for x in fileinput.input(args):
        x = x.rstrip().split("\t")
        seq = x[3]
        w[seq] = w.get(seq, 0) + 1
        l.add(len(seq.split("_")[0]))
    fileinput.close()
    return w, l

def outputMatrix(d, pf, sf, w, lrange, trfs):
    for trf in trfs:
        fstr = "%s.%s.%s.tab" % (pf, trf, sf)
        f = open(fstr, "w")
        trnas = sorted(d[trf].keys())
        f.write("\t".join(["# tRNA"]+map(str, lrange))+"\n")
        for trna in trnas:
            trna_i = { i:0 for i in lrange }
            for seqstr in d[trf][trna]:
                seq, c = seqstr.split("_")
                l, c = len(seq), float(c)/w[seqstr]
                trna_i[l] = trna_i.get(l, 0) + c
            f.write("\t".join([trna]+[str(trna_i[l]) for l in lrange])+"\n")
        f.close()

def main(opts, args):
    if opts.m < 0: raise Exception("bad minimum length")
    if opts.x < 0: raise Exception("bad maximum length")
    if not opts.o:
        pf = os.path.basename(args[0]).split(".")[0]
    else:
        pf = opts.o
    w, l = inputReadStats(args)
    minlen = opts.m if opts.m < sys.maxint else min(l)
    maxlen = opts.x if opts.x > 1 else max(l)
    if maxlen < minlen: raise Exception("error: max_len < min_len")
    lrange = [i for i in xrange(minlen, maxlen+1)]

    trfs = ["5", "3", "l", "t", "i", "o"]
    species = { x:{} for x in trfs }
    isoacceptor = { x:{} for x in trfs }
    isotype = { x:{} for x in trfs }
    for x in fileinput.input(args):
        x = x.rstrip().split("\t")
        trna_i, seqstr, trf = x[0], x[3], x[4].split(":")[0]
        isotype_i, isoacceptor_i = "", ""
        if trna_i.startswith("MT-"):
            trna_i = trna_i.split("_")[0]
            isotype_i = isoacceptor_i = trna_i
        else:
            isotype_i = "-".join(trna_i.split("-")[:-1])
            isoacceptor_i = trna_i.split("_")[0]
        species[trf].setdefault(trna_i, []).append(seqstr)
        isotype[trf].setdefault(isotype_i, []).append(seqstr)
        isoacceptor[trf].setdefault(isoacceptor_i, []).append(seqstr)

    outputMatrix(species, pf, "species", w, lrange, trfs)
    outputMatrix(isotype, pf, "isotype", w, lrange, trfs)
    outputMatrix(isoacceptor, pf, "isoaccp", w, lrange, trfs)


if __name__ == "__main__":
    usage = "%prog [options] tRnaBed"
    description = "Count tRNA fragments by length"

    op = OptionParser(usage=usage, description=description)

    op.add_option("-o", type="string", metavar="PREFIX", default=None,
                  help="Output prefix (default: prefix of input)")
    op.add_option("-m", type="int", metavar="BP", default=sys.maxint,
                  help="Minimum length")
    op.add_option("-x", type="int", metavar="BP", default=1,
                  help="Maximum length")

    (opts, args) = op.parse_args()

    if len(args) == 0: op.error("input BED(s)")

    try:
        main(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
