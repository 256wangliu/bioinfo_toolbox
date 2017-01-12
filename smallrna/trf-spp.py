#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2016 Junko Tsuji

import sys, os.path
from optparse import OptionParser

def main(opts, args):
    nfile, target = args
    lrange, collection = 0, {}
    for x in open(target):
        x = x.rstrip().split("\t")
        if x[0].startswith("#"):
            lrange = len(x[1:])
            print "\t".join(x)
            continue
        collection.update({x[0]:x})
    names = set()
    precursor = ["l", "t", "i"]
    zero = ["0" for i in xrange(lrange)]
    trf, group = target.split(".")[-3:-1]
    for x in open(nfile):
        x, n = x.rstrip(), None
        if x.startswith("MT") or group == "isoaccp":
            n = x.split("_")[0]
        else:
            if group == "species":
                if x.find("_pre_") == -1 and trf not in precursor:
                    n = x
                if x.find("_pre_") >= 0 and trf in precursor:
                    n = x
            elif group == "isotype":
                n = "-".join(x.split("-")[:-1])
        if n: names.add(n)
    for n in sorted(names):
        if n in collection:
            print "\t".join(collection[n])
        else:
            print "\t".join([n] + zero)

if __name__ == "__main__":
    usage = "%prog tRnaNames tRnaLengthMatrix"
    description = "Fill lines for missing tRNAs"

    op = OptionParser(usage=usage, description=description)
    (opts, args) = op.parse_args()

    if len(args) != 2: op.error("input 2 arguments")

    try:
        main(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
