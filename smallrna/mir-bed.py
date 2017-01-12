#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, os.path, fileinput, re
from optparse import OptionParser

def gff3Input(lines):
    Q = [re.compile("Alias="),
         re.compile("Name="),
         re.compile("Derives_from=")]
    for x in lines:
        if x.startswith("#"): continue
        x = x.rstrip().split("\t")
        beg, end, s = int(x[3]), int(x[4]), x[6]
        k = [q.search(x[8]) for q in Q]
        i = x[8][k[0].end():].split(";")[0]
        n = x[8][k[1].end():].split(";")[0]
        p = x[8][k[2].end():].split(";")[0] if k[2] else ""
        yield [x[2], i, beg, end, s, n, p]

def main(opts, args):
    if opts.u < 0 or opts.d < 0:
        raise Exception("input up/down-stream offsets")
    f = {}
    for obj in gff3Input(fileinput.input(args[0])):
        if obj[0] == "miRNA_primary_transcript":
            if obj[4] == "+":
                obj[2] -= opts.u
                obj[3] += opts.d
            else:
                obj[2] -= opts.d
                obj[3] += opts.u
            f.update({obj[1]:obj[2:5]})
    for obj in gff3Input(fileinput.input(args[0])):
        if obj[-1]:
            parent = f[obj[-1]]
            if parent[-1] == "+":
                beg = obj[2] - parent[0]
                end = obj[3] - parent[0] + 1
            else:
                beg = parent[1] - obj[3]
                end = parent[1] - obj[2] + 1
            print "\t".join([obj[-1], str(beg), str(end), obj[-2]])


if __name__ == "__main__":
    usage = "%prog [options] gff3File"
    description = "Extract mature miRNA annotation and convert the genomic loci to relative positions based on miRNA primary transcripts"

    op = OptionParser(usage=usage, description=description)
    op.add_option("-u", type="int", metavar="BP", default=-1,
                  help="Upstream offset")
    op.add_option("-d", type="int", metavar="BP", default=-1,
                  help="Downstream offset")

    (opts, args) = op.parse_args()

    if len(args) != 1: op.error("input a gff3 file")

    try:
        main(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
