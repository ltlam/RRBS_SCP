#=======================================================================================================================
# bam_to_cov.py - Script to report the coverage for mappable regions
# Created By: Larry Lam
# Date Created: 12/16/2015
#=======================================================================================================================
from optparse import OptionParser
import pysam
import numpy as np

STR_ALL = 'all'

#-----------------------------------------------------------------------------------------------------------------------
# getMapRegions
# [in] read map regions file for RRBS
#-----------------------------------------------------------------------------------------------------------------------
def getMapRegions(strMapRegFile,strChrSel):
    lfrags = []
    for strLine in open(strMapRegFile,'r'):

        if strMapRegFile.endswith('.txt'):
            strChr, strRow, strStart, strEnd, strFrag = strLine.rstrip().split()
        else:
            strChr, strStart, strEnd = strLine.rstrip().split()

        if strChr != strChrSel and strChrSel != STR_ALL:
            continue
        lfrags.append((strChr,int(strStart),int(strEnd)))
    return lfrags

#-----------------------------------------------------------------------------------------------------------------------
# rptCover
# [in] dFragInfo - fragments
# [in] strBamFile - path to bam file
#-----------------------------------------------------------------------------------------------------------------------
def rptCover(lfrags,strBamFile,strChrSel):
    if strChrSel != STR_ALL:
        str_out_fname = strBamFile+'.%s.cov'%strChrSel
    else:
        str_out_fname = strBamFile+'.cov'

    fout = open(str_out_fname,'w')
    samfile = pysam.Samfile( strBamFile, 'rb')
    n_cover_sum = 0
    n_pos_count = 0
    strCurChr = ''
    for i in range(len(lfrags)):
        strChr,nStart,nEnd = lfrags[i]
        n_pos_count += (nEnd - nStart + 1)
        if strCurChr != strChr:
            print strChr
            strCurChr = strChr
        for pileupcolumn in samfile.pileup(strChr, nStart, nEnd):
            n_cover_sum += pileupcolumn.n
    samfile.close()
    fout.write('%d\t%d\t%f\n'%(n_cover_sum,n_pos_count,float(n_cover_sum)/float(n_pos_count)))
    fout.close()

#-----------------------------------------------------------------------------------------------------------------------
# main
#-----------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-r", "--reg", dest="reg", help="RRBS mappable region file", metavar="REGF")
    parser.add_option("-b", "--bam", dest="bam", help="bam file", metavar="BAM")
    parser.add_option("-c", "--chr", dest="chr", help="chr chromosome", metavar="CHR",default=STR_ALL)

    (options, args) = parser.parse_args()

    #getMapRegion Info
    lfrags = getMapRegions(options.reg,options.chr)

    #printCoverage
    rptCover(lfrags,options.bam,options.chr)

    print '--END--'