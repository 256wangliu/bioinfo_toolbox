#! /usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2015 Junko Tsuji

# Calculates molecular weights of amino acid sequences
# from FASTA files. Molecular weights are computed with
# the following formula:
#     X = Sum(Weights of individual residues)
#     H = Water weight (18.015)
#     C = Number of regidues - 1
#     Molecular weight = X - H * C

import sys, os.path, fileinput
from optparse import OptionParser

H2O = 18.015
aa_mw = {'A': 89.094, 'B':132.611, 'C':121.159, 'D':133.103,
         'E':147.130, 'F':165.190, 'G': 75.067, 'H':155.155,
         'I':131.174, 'K':146.188, 'L':131.174, 'M':149.212,
         'N':132.118, 'P':115.131, 'Q':146.145, 'R':174.202,
         'S':105.093, 'T':119.120, 'V':117.147, 'W':204.226,
         'Y':181.189, 'Z':146.638, 'X':128.160}

def fastaInput(lines):
    name, seq = "", ""
    for line in lines:
        line = line.rstrip()
        if not line: continue
        if line.startswith('>'):
            if seq != "": yield name, seq
            name = line.replace('>', '')
            seq = ""
        else:
            seq += line.upper()
    if seq != "": yield name, seq

def computeMolecularWeight(seq):
    seq = seq.replace('*', '')
    seqLen = len(seq)
    sumWeight = sum([aa_mw[c] for c in seq])
    return (sumWeight - (H2O * (seqLen - 1)))/1000

def main(opts, args):
    fasta = fastaInput(fileinput.input(args))
    for name, seq in fasta:
        mw = computeMolecularWeight(seq)
        if opts.l <= mw and mw <= opts.u:
            print '\t'.join([name, str(round(mw, 2))])

if __name__ == "__main__":
    usage = "%prog fastaFile(s)"
    description = "Compute molecular weight (kDa) from amino acid sequences in FASTA"

    op = OptionParser(usage=usage, description=description)

    op.add_option("-l", metavar="kDa", type=float, default=0,
                  help="Lower bound of kDa")
    op.add_option("-u", metavar="kDa", type=float, default=sys.maxint,
                  help="Upper bound of kDa")

    (opts, args) = op.parse_args()

    if len(args) < 1: op.error("input FASTA file(s)")

    try:
        main(opts, args)
    except KeyboardInterrupt: pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
