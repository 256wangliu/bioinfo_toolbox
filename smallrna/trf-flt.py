#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2016 Junko Tsuji

import sys, os.path, fileinput
from optparse import OptionParser

def getInfo(seq, off_5, off_3, his, annot, overlap):
    L = len(seq.split("_")[0])
    if annot == "mature_trna":
        if off_5 <= 5 + his:
            nucleotide = ":G" if off_5 == 0 and his else ""
            return ("5:%d" + nucleotide) % L, str(off_5)
        elif -7 <= off_3:
            nucleotide = ":CCA"[::-1][abs(off_3):][::-1] if -2 <= off_3 else ""
            return ("3:%d" + nucleotide) % L, str(off_3)
        else:
            return "o:%d" % L, "."
    elif annot == "5_leader" and abs(off_3) <= 5:
        return "l:%d" % L, str(off_3)
    elif annot == "3_trailer" and abs(off_5) <= 5:
        return "t:%d" % L, str(off_5)
    elif annot == "intron" and int(overlap) >= 5:
        return "i:%d" % L, overlap
    else:
        return None, None

def main(opts, args):
    stats = {}
    collection = {}
    isPseudo = lambda o: o.startswith("Pseudo")
    isPrecursor = lambda o: o.find("_pre_") >= 0
    for x in fileinput.input(args):
        x = x.rstrip().split("\t")
        seq = x[3]
        his = x[0].startswith("His")
        off_5 = int(x[1])-int(x[7])
        off_3 = int(x[2])-int(x[8])
        flag, nt = getInfo(seq, off_5, off_3, his, x[9], x[10])
        if not flag and not nt: continue
        collection.setdefault(seq, []).append(x[:4]+[flag, nt])
        i = 0
        if isPrecursor(x[0]): i = 1
        elif isPseudo(x[0]):  i = 2
        stats.setdefault(seq, [0,0,0])[i] += 1
    for seq in collection:
        func = None
        if stats[seq][0] > 0: # print only mature tRNAs
            func = lambda o: not isPseudo(o) and not isPrecursor(o)
        else:
            if stats[seq][1] > 0: # exclude pseudogenes
                func = lambda o: not isPseudo(o)
            else:
                func = lambda o: True
        for x in collection[seq]:
            if func(x[0]):
                print "\t".join(x)

if __name__ == "__main__":
    usage = "%prog [options] intersectedBed(s)"
    description = "Filter, assign, and annotate tRNA fragments"

    op = OptionParser(usage=usage, description=description)
    (opts, args) = op.parse_args()

    if len(args) == 0: op.error("input BED(s)")

    try:
        main(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
