#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2016 Junko Tsuji

import sys, os.path, fileinput, string
from optparse import OptionParser

codon = { "TTT":"Phe", "TTC":"Phe", "TTA":"Leu", "TTG":"Leu",
          "CTT":"Leu", "CTC":"Leu", "CTA":"Leu", "CTG":"Leu",
          "ATT":"Ile", "ATC":"Ile", "ATA":"Ile", "ATG":"Met",
          "GTT":"Val", "GTC":"Val", "GTA":"Val", "GTG":"Val",
          "TCT":"Ser", "TCC":"Ser", "TCA":"Ser", "TCG":"Ser",
          "CCT":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
          "ACT":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",
          "GCT":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",
          "TAT":"Tyr", "TAC":"Thr", "TAA":"Sup", "TAG":"Sup",
          "CAT":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",
          "AAT":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",
          "GAT":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
          "TGT":"Cys", "TGC":"Cys", "TGA":"SeC", "TGG":"Trp",
          "CGT":"Arg", "CGC":"Arg", "CGA":"Arg", "CGG":"Arg",
          "AGT":"Ser", "AGC":"Ser", "AGA":"Arg", "AGG":"Arg",
          "GGT":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly"  }

complement = string.maketrans("ACGTRYKMBDHVacgtrykmbdhv",
                              "TGCAYRMKVHDBtgcayrmkvhdb")
def reverseComplement(seq):
    return seq[::-1].translate(complement)

def main(args):
    for x in fileinput.input(args[0]):
        x = x.rstrip("\n").replace(" ","").split("\t")
        chrom, beg, end = x[0], int(x[2]), int(x[3])
        if beg > end:
            strand, beg, end = "-", end, beg
        else:
            strand = "+"
        beg -= 1
        aa = x[4]
        anticodon = x[5].upper()
        if aa == "Pseudo":
            flag = "1"
            if anticodon != "???":
                aa = aa + "-" + codon[reverseComplement(anticodon)]
        else:
            flag = "0"
        name = aa + "-" + anticodon
        intron_beg, intron_end = int(x[6]), int(x[7])
        if (intron_beg + intron_end) > 0:
            if strand == "-":
                intron_beg, intron_end = intron_end, intron_beg
            intron_beg -= 1
            intron = "intron:" + "-".join(map(str, [intron_beg, intron_end]))
        else:
            intron = "no_intron"
        print "\t".join([chrom, str(beg), str(end), name, flag, strand, intron])


if __name__ == "__main__":
    usage = "%prog tRNAscanOutput"
    description = "Convert tRNAscan-SE position output to BED"

    op = OptionParser(usage=usage, description=description)
    (opts, args) = op.parse_args()

    if len(args) != 1: op.error("input tRNAscan-SE output")

    try:
        main(args)
    except KeyboardInterrupt: pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
