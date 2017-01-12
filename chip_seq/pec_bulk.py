#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2016 Junko Tsuji

# To get the Synapse IDs of datasets, the following commands were run:
#
#   $ synapse query 'SELECT id,disease,tissueType,assayTarget,cellType FROM file \
#                    WHERE projectId=="syn4921369" and dataType=="histone modification" and \
#                          organism=="Homo sapiens" and fileType=="bam" and study=="EpiMap"' | \
#     grep Prefrontal | grep NeuN | grep H3K4me3 | cut -f1 > ID_LIST

import sys, os.path, fileinput, subprocess
from optparse import OptionParser
import synapseclient

def main(args):
    # beforehand, make ~/.synapseConfigure
    syn = synapseclient.Synapse()
    syn.login()
    path = args[0]
    for x in fileinput.input(args[1:]):
        xid = x.rstrip().split()[0]
        info = syn.get(xid, downloadFile=False)
        F = path + "/" + info.name
        if os.path.exists(F):
            sys.stderr.write("%s:%s exists, skipped\n" % (xid, F))
            continue
        obj = syn.get(xid, downloadLocation=path)
        md5_xid = subprocess.check_output(["md5sum", F]).split()[0]
        if md5_xid != obj.md5:
            sys.stderr.write("%s:%s md5sum failed\n" % (xid, F))
            subprocess.call(["rm", F])
        else:
            sys.stderr.write("%s:%s successfully downloaded\n" % (xid, F))
    syn.logout()

if __name__ == "__main__":
    usage = "%prog downloadPath idList(s)"
    description = "Download BAM files at once by specifying Synapse IDs"

    op = OptionParser(usage=usage, description=description)
    (opts, args) = op.parse_args()

    try:
        main(args)
    except KeyboardInterrupt: pass
    except Exception, e:
        prog = os.path.basename(sys.argv[0])
        sys.exit(prog + ": error: " + str(e))
