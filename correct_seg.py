#!/usr/bin/env python
# author: Junko Tsuji <jtsuji@broadinstitute.org>
# this script is written based on Justin Rhoades' python script


import sys
import math
import os.path
from argparse import ArgumentParser


def main(args):
    fo = open(args.output, "w") if args.output else sys.stdout

    # diploid
    CN = 2

    TF = args.purity
    PL = args.ploidy

    copy_ratio_index = None
    copy_status_index = None

    for i, x in enumerate(open(args.input_seg)):
        x = x.rstrip().split("\t")
        if i == 0:
            copy_ratio_index = x.index(args.copy_ratio_column)
            if args.copy_status_column:
                copy_status_index = x.index(args.copy_status_column)
            fo.write("\t".join(x) + "\n")
            continue

        # normalize tumor copy ratio with purity and ploidy
        if x[copy_ratio_index] != "NA":
            ct = (2**float(x[copy_ratio_index]) * (CN*(1-TF) + TF*PL) - CN*(1-TF)) / TF
            normalized_ratio = max(ct, 1.0/(2**6)) / PL
            x[copy_ratio_index] = str(math.log(normalized_ratio, 2))

        # if copy status column is specified, set neutral segment copy ratios to zero
        if copy_status_index != None:
            x[copy_ratio_index] = "0.0" if x[copy_status_index] == "NEUT" else x[copy_ratio_index]

        fo.write("\t".join(x) + "\n")
    fo.close()


if __name__ == "__main__":
    parser = ArgumentParser(description="Correct purity and ploidy of a sample")

    parser.add_argument("--purity", type=float, help="tumor fraction", required=True)
    parser.add_argument("--ploidy", type=float, help="sample ploidy", required=True)
    parser.add_argument("--copy-ratio-column", type=str,
                        help="name of copy ratio column, e.g. 'median', 'seg.median.log2R'", required=True)

    parser.add_argument("--copy-status-column", type=str,
                        help="name of copy status column to set neutral segment copy ratios to zero, e.g. 'event', 'call' (default: deactivated)", default=None)
    parser.add_argument("--output", help="output file name (default: stdout)", default=None)

    parser.add_argument("input_seg", help="input seg file")
    
    args = parser.parse_args()
    try:
        main(args)
    except KeyboardInterrupt: pass

