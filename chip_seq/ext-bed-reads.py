#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os.path, re, fileinput
from optparse import OptionParser

def main(opts, args):
    extLen = 0
    p = re.compile("[MIDNSHP=X]")
    try:
        extLen = int(args[1])
    except TypeError:
        raise Exception("bad extension length")
    if extLen <= 0:
        raise Exception("bad extension length")
    if not os.path.exists(args[0]):
        raise Exception("can't read " + args[0])
    if not os.path.exists(args[2]):
        raise Exception("can't read " + args[2])
    chrom = {}
    for c in fileinput.input(args[0]):
        c = c.rstrip().split()
        chrom.update({c[0]:int(c[1])})
    prefix = ("EXTEND_FRAG=%dbp_" % extLen)
    for x in fileinput.input(args[2]):
        x = x.rstrip().split()
        xExtLen = extLen
        cigar = p.findall(x[-1])
        cilen = map(int, p.split(x[-1])[:-1])
        for i in xrange(len(cigar)):
            if   cigar[i] == 'M': xExtLen -= cilen[i]
            elif cigar[i] == 'I': xExtLen -= cilen[i]
            elif cigar[i] == '=': xExtLen -= cilen[i]
            elif cigar[i] == 'X': xExtLen -= cilen[i]
        if xExtLen < 0:
            raise Exception("specify larger extension length")
        if x[5] == '-':
            x[1] = str(max(int(x[1])-xExtLen, 0))
        else:
            x[2] = str(min(int(x[2])+xExtLen, chrom[x[0]]))
        x[3], x[4] = (prefix + x[3]), '.'
        print "\t".join(x[:-1])

if __name__ == "__main__":
    usage = "%prog chromLen extLen mappedReadBed"
    description = "Extend mapped ChIP-seq reads"

    op = OptionParser(usage=usage, description=description)
    (opts, args) = op.parse_args()

    if len(args) != 3: op.error("input three arguments")

    try:
        main(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
