#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os.path, fileinput
from optparse import OptionParser

def filter(group, i):
    m = min([abs(int(x[i])) for x in group])
    return [x for x in group if m == abs(int(x[i]))]

def assignMirRead(group):
    if len(group) == 1:
        print "\t".join(group[0])
    else:
        g_opt = filter(group, -2)
        g_opt = filter(g_opt, -1)
        for x in g_opt:
            print "\t".join(x)

def main(opts, args):
    sp = ""
    group = []
    for x in fileinput.input(args):
        x = x.rstrip().split("\t")
        if sp != x[3]:
            if not sp: pass
            else:      assignMirRead(group)
            sp = x[3]
            group = [x]
        else:
            group.append(x)
    assignMirRead(group)


if __name__ == "__main__":
    usage = "%prog sortedMirBed"
    description = "Assign multimply mapped reads to miRNAs that have minimum offsets"

    op = OptionParser(usage=usage, description=description)
    (opts, args) = op.parse_args()

    if len(args) != 1: op.error("input a bed input")

    try:
       main(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
         prog = os.path.basename(sys.argv[0])
         sys.exit(prog + ": error: " + str(e))
