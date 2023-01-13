#!/usr/bin/env python
"""Script to collect and calcuate the frips statistics from across all runs. 
#Total Peaks, #10FC peaks, #20FC peaks
Outputs: 
Run,Total,ReadsInPeaks,FRiP
"""

import os
import sys
from optparse import OptionParser

def main():
    usage = "USAGE: %prog -f [FILE_1] -f [FILE_2] ...-f [FILE_N] -o [output]"
    optparser = OptionParser(usage=usage)
    optparser.add_option("-f", "--files", action="append", help="list of _sorted_peaks.narrowPeak files")
    optparser.add_option("-o", "--output", help="output file")

    (options, args) = optparser.parse_args(sys.argv)

    if not options.files or not options.output:
        optparser.print_help()
        sys.exit(-1)

    out = open(options.output,"w")
    out.write("%s\n" % ",".join(["Run","Total","ReadsInPeaks","FRiP"]))

    for f in options.files:
        #TRY to infer the RUN NAMES
        runID = ".".join(f.strip().split("/")[-1].split('.')[:-1])
        if runID.endswith('_frip'):
            runID = runID.replace("_frip","")

        f = open(f)
        readsInPeaks = int(f.readline().strip().split("\t")[1])
        totalPeaks = int(f.readline().strip().split("\t")[1])
        frip = "%.1f" % (float(readsInPeaks)/totalPeaks *100.0)

        out.write("%s\n" % ",".join([runID,str(totalPeaks),str(readsInPeaks),str(frip)]))
        f.close()
    out.close()

if __name__=='__main__':
    main()


