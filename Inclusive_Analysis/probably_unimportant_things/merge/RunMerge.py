#!/usr/bin/python

from glob import glob
from subprocess import Popen, PIPE
import os

#_____________________________________________________________________________
if __name__ == "__main__":

    #merge output files by chunks of a given size

    top = "/gpfs01/star/pwg/jaroslav/test/star-upcDst/trees/UPC_main_JpsiB_10_11_14_v2"
    pattern = "/{r14_prod,r14_low,r14_mid,r14_high}/*.root"
    #pattern = "/r11/*.root"

    outdir = "merge_r14"
    outfile = "StUPC_main_JpsiB_14.root"
    outlist = "StUPC_main_JpsiB_14.list"

    chunksiz = int(2.5e6)  # approx kB

    #-- end of config --


    #create output directory for merged files
    merge_path = top+"/"+outdir
    if os.access(merge_path, os.F_OK) == False:
        os.makedirs(merge_path)

    #output files list
    chunk_list = open(merge_path+"/"+outlist, "w")

    #list done jobs
    totsiz = 0    # total chunk size
    ichunk = 0    # chunk index
    nchunk = 0    # number of files in chunk
    tmpnam = "files.tmp"   # temporary to list files for a given chunk
    tmp = open(tmpnam, "w")
    cmd = "ls -s " + top + pattern
    out = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()[0].split("\n")
    print "Number of all outputs:", len(out)
    print
    #list loop
    for iline in xrange(len(out)):
        fline = out[iline]
        size = 0
        if len(fline) > 0:
            siznam = fline.lstrip().split(" ")
            size = int(siznam[0])
            name = siznam[1]
        if size > 0:
            totsiz += size
            tmp.write(name+"\n")
            nchunk += 1
        if totsiz >= chunksiz or iline == len(out)-1:
            tmp.close()
            print "Chunk:", ichunk
            print "Num files:", nchunk
            print "Size:", totsiz
            iform = "_{0:04d}".format(ichunk)
            chunkout = merge_path+"/"+outfile.split(".root")[0]+iform+".root"
            merge_cmd = "root -l -b -q MergeFiles.C(\"" + tmpnam + "\",\"" + chunkout + "\")"
            merg = Popen(merge_cmd.split(), stdout=PIPE, stderr=PIPE).communicate()
            print merg[0], merg[1]
            chunk_list.write(chunkout+"\n")
            #reset the temporary and indices
            tmp = open(tmpnam, "w")
            totsiz = 0
            nchunk = 0
            #increment the chunk index
            ichunk += 1
    #list loop end

    #clean-up the list temporary file
    tmp.close()
    os.remove(tmpnam)

    print "All done."
























