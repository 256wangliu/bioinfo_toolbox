#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2015 Junko Tsuji

# This generates a BED file of precursor transcripts from a GFF3
# annotation provided by miRBase: http://www.mirbase.org/ftp.shtml
# To allow isomir search, up/down- stream values can be accepted.

import sys, os.path, fileinput
from optparse import OptionParser

def main(opts, args):
    if opts.u < 0: raise Exception("bad value: -u")
    if opts.d < 0: raise Exception("bad value: -d")
    for x in fileinput.input(args[0]):
        x = x.rstrip()
        if x.startswith("#") or not x: continue
        x = x.split("\t")
        if x[2] != "miRNA_primary_transcript":
            continue
        chrom = x[0]
        beg = int(x[3]) - 1
        end = int(x[4])
        s = x[6]
        name = x[8].split("Alias=")[1].split(";")[0]
        if s == "-":
            beg -= opts.d
            end += opts.u
        else:
            beg -= opts.u
            end += opts.d
        print "\t".join([chrom, str(beg), str(end), name, ".", s])

if __name__ == "__main__":
    usage = "%prog [option] gff3File"
    description = "Generate BED of miRNA primary transcripts from GFF3 provided by mirBase"
    
    op = OptionParser(usage=usage, description=description)
    op.add_option("-u", type="int", default=5, metavar="BP",
                  help="Upstream length to include (default=%default)")
    op.add_option("-d", type="int", default=5, metavar="BP",
                  help="Downstream length to include (default=%default)")

    (opts, args) = op.parse_args()

    if len(args) != 1: op.error("input gff3 file")
    try:
        main(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
