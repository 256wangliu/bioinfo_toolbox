#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os.path, fileinput
from optparse import OptionParser

def main(opts, args):
    if opts.x < 0: raise Exception("bad maximum length")
    if opts.m < 0: raise Exception("bad minimum length")
    d, s = {}, {}
    for x in fileinput.input(args[0]):
        x = x.rstrip().split()
        l, c = len(x[0]), float(x[1])
        d[l] = d.get(l, 0) + c
        s[l] = s.get(l, 0) + 1
    L = d.keys()
    minlen = opts.m if opts.m < sys.maxint else min(L)
    maxlen = opts.x if opts.x > 1 else max(L)
    # output will be: (1) length  (2) nonredundant_counts  (3) read_counts
    for i in xrange(minlen-1, maxlen):
        l = i + 1
        if l in d:
            print "\t".join([str(l), str(s[l]), str(d[l])])
        else:
            print "\t".join([str(l), "0", "0.0"])

if __name__ == "__main__":
    usage = "%prog [options] readCounts"
    description = "Output read length distribution. Input has to be composed two columns: (1) read sequence and (2) its count"

    op = OptionParser(usage=usage, description=description)
    op.add_option("-m", type="int", metavar="BP", default=sys.maxint,
                  help="Minimum length")
    op.add_option("-x", type="int", metavar="BP", default=1,
                  help="Maximum length")

    (opts, args) = op.parse_args()

    if len(args) != 1: op.error("input 1 argument")

    try:
        main(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
