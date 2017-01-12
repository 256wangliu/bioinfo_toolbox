#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os.path, fileinput
from optparse import OptionParser

def main(opts, args):
    if opts.m < 0: raise Exception("bad minimum length")
    if opts.x < 0: raise Exception("bad maximum length")
    if opts.m >= opts.x: raise Exception("bad length range")
    collection = {}
    trfs = ["5", "3", "l", "t", "i", "o"]
    index = { x:i for i,x in enumerate(trfs) }
    for x in fileinput.input(args):
        x = x.rstrip().split("\t")
        trf = x[4].split(":")[0]
        collection.setdefault(x[3], [0,0,0,0,0,0])[index[trf]] += 1
    a = {}  # all reads
    d = { i:{} for i in trfs }
    for seqstr in collection:
        seq, c = seqstr.split("_")
        l = len(seq)
        counts = collection[seqstr]
        unit = float(c) / sum(counts)
        for i, trf in enumerate(trfs):
            if counts[i] == 0: continue
            share = unit * counts[i]
            d[trf][l] = d[trf].get(l, 0) + share
        a[l] = a.get(l, 0) + int(c)
    obs_len = a.keys()
    minlen = opts.m if opts.m < sys.maxint else min(obs_len)
    maxlen = opts.x if opts.x > 1 else max(obs_len)
    if maxlen < minlen: raise Exception("error: max_len < min_len")

    header = ["# length", "total"] + trfs
    print "\t".join(header)
    for l in xrange(minlen, maxlen+1):
        range_i = [str(d[trf].get(l, 0)) for trf in trfs]
        print "\t".join([str(l), str(a.get(l, 0))] + range_i)

if __name__ == "__main__":
    usage = "%prog [optinos] tRnaBed"
    description = "Read length distribution of tRNA fragments"

    op = OptionParser(usage=usage, description=description)

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
