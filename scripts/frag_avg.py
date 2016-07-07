import argparse
import numpy as np
import re

N_SAM_COL = 2

def load_fragmeth(str_fragmeth):
    d_frag_to_col_to_lmeth = {}
    l_header = []
    for str_line in open(str_fragmeth):
        lcols = str_line.rstrip().split('\t')
        if len(l_header) == 0:
            l_header = lcols
            continue
        str_frag = lcols[1]
        n_sam = len(lcols) - N_SAM_COL
        for i in range(N_SAM_COL,len(lcols)):
            if d_frag_to_col_to_lmeth.has_key(str_frag) == False:
                d_frag_to_col_to_lmeth[str_frag] = {}
                for j in range(N_SAM_COL,len(lcols)):
                    d_frag_to_col_to_lmeth[str_frag][j] = []
            d_frag_to_col_to_lmeth[str_frag][i].append(float(lcols[i]))
    return(l_header,d_frag_to_col_to_lmeth)

def write_frag_avg(lheader,d_frag_to_col_to_lmeth,str_out_path):
    fout = open(str_out_path,'w')

    str_header = 'FRAG\tNUM_CG'
    for i in range(N_SAM_COL,len(lheader)):
        str_header+='\t%s'%lheader[i]
    fout.write(str_header+'\n')

    for str_frag in sorted(d_frag_to_col_to_lmeth.keys()):
        str_out_line = '%s\t%d'%(str_frag,len(d_frag_to_col_to_lmeth[str_frag][N_SAM_COL]))
        for i in range(N_SAM_COL,len(lheader)):
            str_out_line+='\t%0.12f'%(float(sum(d_frag_to_col_to_lmeth[str_frag][i]))/float(len(d_frag_to_col_to_lmeth[str_frag][i])))
        fout.write(str_out_line+'\n')
    fout.close()

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--ffcgmap', type=str, help="sample pos and fragment methylation")
    args = parser.parse_args()

    #read ffcgmap
    lheader, d_frag_to_col_to_lmeth = load_fragmeth(args.ffcgmap)

    #write frag_avg
    write_frag_avg(lheader,d_frag_to_col_to_lmeth,re.sub('_frag.txt$','_frag_avg.txt',args.ffcgmap))

