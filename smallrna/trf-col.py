#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2016 Junko Tsuji

import sys, os.path
from optparse import OptionParser, SUPPRESS_HELP

def main(opts, args):
    if opts.m < 0: raise Exception("bad minimum length")
    if opts.x < 0: raise Exception("bad maximum length")
    if opts.m >= opts.x: raise Exception("bad length range")

    ml, xl = opts.m, opts.x
    ml_i, xl_i = -1, -1

    rownames = []
    for x in open(args[0]):
        x = x.rstrip().split("\t")
        if x[0].startswith("# tRNA"):
            l = map(int, x[1:])
            ml_i = l.index(ml)
            xl_i = l.index(xl)
            continue
        rownames.append(x[0])

    mat, samples = [], []
    for f in args:
        sample = os.path.basename(f).split(opts.d)[0]
        samples.append(sample)
        l = [ sum(map(float, x.rstrip().split("\t")[1:])[ml_i : xl_i+1])
              for x in open(f) if not x.startswith("#") ]
        mat.append(map(str, l))

    print "\t".join(["# tRNA"] + samples)
    for trna in xrange(len(rownames)):
        l = [ mat[i][trna] for i in xrange(len(samples)) ]
        print "\t".join([rownames[trna]] + l)


if __name__ == "__main__":
    usage = "%prog [options] LengthMatrices"
    description = "Make read count table in a specific length range"

    op = OptionParser(usage=usage, description=description)

    op.add_option("-m", type="int", metavar="BP", default=sys.maxint,
                  help="Minimum length to extract")
    op.add_option("-x", type="int", metavar="BP", default=1,
                  help="Maximum length to extract")
    op.add_option("-d", type="string", metavar="DELIM", default=".",
                  help=SUPPRESS_HELP)
    (opts, args) = op.parse_args()

    try:
        main(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
