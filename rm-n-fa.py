#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2015 Junko Tsuji

# This simply removes FASTA sequences composed of only 'N'.

import sys, os.path, fileinput
from optparse import OptionParser

def fastaInput(lines):
    tag, seq = "", ""
    for x in lines:
        x = x.rstrip()
        if x.startswith(">"):
            if seq: yield tag, seq
            tag, seq = x, ""
        else:
            seq += x.upper()
    if seq: yield tag, seq

def main(opts, args):
    for tag, seq in fastaInput(fileinput.input(args[0])):
        if len(seq.replace('N','')) == 0: continue
        print '\n'.join([tag, seq])

if __name__ == "__main__":
    usage = "%prog [option] fastaFile"
    description = "Remove sequences composed of 'N' bases from FASTA"
    
    op = OptionParser(usage=usage, description=description)
    (opts, args) = op.parse_args()

    if len(args) != 1: op.error("input FASTA")
    try:
        main(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
