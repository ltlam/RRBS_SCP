#==============================================================================
# cgmap_data_frame.py
# Created By: Larry Lam
# Date Created: 3/17/2016
# 1. sort bed
# 2. perform union
#==============================================================================
from multiprocessing import Pool
import pandas as pd
import numpy as np
import argparse
import shutil
import re
import os

STR_TMP_DIR = 'decon_tmp_cgmap'
STR_TMP_CSV = 'merge_cgmap.txt'
STR_CNTXT_SEL = 'CG'
N_MIN_COV = 10
STR_MAX_CHR = 'chr9999'

#------------------------------------------------------------------------------
# sort_cgmap_wrapper - wrapper for sort select cgmap
# [in] args - sort_cgmap wrapper
#------------------------------------------------------------------------------
def sort_cgmap_wrapper(args):
    return sort_cgmap(*args)

#------------------------------------------------------------------------------
# sort_cgmap - sort cgmap to ascending by chr and position and
#                 write to temp folder
# [in] str_cgmap - path to bed graph file
# [in] str_temp_cgmap_folder - path to temporary folder of filtered cgmap
# [in] str_cntxt - context filter
# [in] n_min_cov - minimum coverage
# [in] n_sorted - flag 1 = sort, 0 = pre sorted
#------------------------------------------------------------------------------
def sort_cgmap(str_cgmap,str_temp_cgmap_folder,str_cntxt,n_min_cov,n_sorted):
    #check if temp folder exists
    head,tail = os.path.split(str_cgmap)
    str_tmp_fld = os.path.join(head,str_temp_cgmap_folder)
    str_sorted_file = os.path.join(str_tmp_fld,tail)
    if not os.path.exists(str_tmp_fld):
        os.makedirs(str_tmp_fld)

    #if sorted copy to temp folder
    if not n_sorted:
        cgmap_df = pd.read_table(str_cgmap,sep='\t',header=None)
        sorted_df = cgmap_df.sort_values(by=[0,2],ascending=[1,1])
        if str_cntxt != '':
            sorted_df = sorted_df.loc[sorted_df.iloc[:,3]==str_cntxt]
        sorted_df = sorted_df[sorted_df.iloc[:,7] >= n_min_cov]
        sorted_df.to_csv(str_sorted_file,na_rep='NA',sep='\t',header=False,
                         index=False)
    else:
        shutil.copyfile(str_cgmap, str_sorted_file)

#------------------------------------------------------------------------------
# write_intersection - merge the methylation calls into a file based data frame
# [in] l_sorted_cgmap - list of the cgmap sorted by pos
# [in] str_out - path to output file
#------------------------------------------------------------------------------
def write_intersection(l_sorted_cgmap,str_out):
    #init temp variables
    l_cur_chr = []
    l_cur_pos = []
    l_cur_meth = []
    fout = open(str_out,'w')

    #open cgmap files
    l_fopen = []
    for str_cur_bed in l_sorted_cgmap:
        l_fopen.append(open(str_cur_bed,'r'))

    #init cur pos
    if len(l_cur_chr) == 0:
        for n_fid in range(len(l_fopen)):
            str_file_line = l_fopen[n_fid].readline()
            lcols = str_file_line.split('\t')
            l_cur_chr.append(lcols[0])
            l_cur_pos.append(int(lcols[2]))
            l_cur_meth.append(float(lcols[5]))

    #write header
    str_temp_row = 'POS'
    for i in range(len(l_sorted_cgmap)):
        head,tail =os.path.split(l_sorted_cgmap[i])
        str_temp_row+='\t%s'%re.sub('\\.txt$|\\.bed$|\\.CGmap$','',tail)
    fout.write(str_temp_row+'\n')

    while(l_cur_chr.count(STR_MAX_CHR) != len(l_sorted_cgmap)):
        #init meth val
        l_temp_meth = [float('nan')] * len(l_sorted_cgmap)

        #check min id
        str_min_chr = min(l_cur_chr)
        l_min_chr_id = [i for i, x in enumerate(l_cur_chr) if x==str_min_chr]
        l_temp_pos = l_cur_pos[:]
        for i in range(len(l_temp_pos)):
            if not (i in l_min_chr_id):
                l_temp_pos[i] = float('inf')

        #check min pos
        n_min_pos = min(l_temp_pos)
        l_min_pos_id = [i for i, x in enumerate(l_temp_pos) if x==n_min_pos]
        if n_min_pos==float('inf'):
            break

        #assign meth to temp list
        for n_cur_min_id in l_min_pos_id:
            l_temp_meth[n_cur_min_id] = l_cur_meth[n_cur_min_id]

        str_temp_row = '%s_%d'%(str_min_chr,n_min_pos)
        for i in range(len(l_temp_meth)):
            str_temp_row+= '\t' + '%.12f'%l_temp_meth[i]

        #if temp meth contains nan, continue
        if not np.isnan(l_temp_meth).any():
            fout.write(str_temp_row+'\n')

        #update chr, pos, curmeth, and fin
        for n_min_id in l_min_pos_id:
            str_file_line = l_fopen[n_min_id].readline()
            if str_file_line != '':
                lcols = str_file_line.split('\t')
                l_cur_chr[n_min_id] = lcols[0]
                l_cur_pos[n_min_id] = int(lcols[2])
                l_cur_meth[n_min_id] = float(lcols[5])
            else:
                l_cur_chr[n_min_id] = STR_MAX_CHR
                l_cur_pos[n_min_id] = float('inf')
                l_cur_meth[n_min_id] = float('nan')

    #close bedgraph files
    for fopen in l_fopen:
        fopen.close()
    fout.close()

def get_parser(parser=None):
    if parser==None:
        parser = argparse.ArgumentParser()
    parser.add_argument('-in_sam', '--in_sam', type=str,
                        help="comma separated list of cgmap files")
    parser.add_argument('-out_sam','--out_sam',
                        type=str,help="output file name",
                        default=STR_TMP_CSV)
    parser.add_argument('-f_sam', '--filter_context', type=str,
                        help="filter context CG,CHG,CHH", default=STR_CNTXT_SEL)
    parser.add_argument('-cov_sam', '--min_cov', type=int,
                        help="minimum cover threshold", default=N_MIN_COV)
    parser.add_argument('-s_sam', '--sorted', type=int,
                        help="sort file 0=unsorted, 1=sorted", default=0)
    parser.add_argument('-p_sam','--proc_sam',
                        type=int,
                        help="number of threads (default 1)",
                        default=1)
    return parser

def get_files_to_sam_lbl(str_cgmap_files):
    l_sams = []
    l_cgmap_files = str_cgmap_files.split(',')
    for str_cgmap in l_cgmap_files:
        head,tail = os.path.split(str_cgmap)
        l_sams.append(re.sub('.CGmap$','',tail))
    return ','.join(l_sams)

def main(str_cgmaps,n_min_cov,str_out,str_context,n_sorted,n_process=1):
    #sort files
    l_sorted_cgmap = []
    l_multi_args = []
    l_cgmap_files = str_cgmaps.split(',')
    for str_cgmap in l_cgmap_files:
        head,tail = os.path.split(str_cgmap)
        l_sorted_cgmap.append(os.path.join(head,STR_TMP_DIR+str_out.split('.')[0],tail))
        l_multi_args.append([str_cgmap,STR_TMP_DIR+str_out.split('.')[0],str_context,n_min_cov,n_sorted])

    '''
    pool = Pool(n_process)
    pool.map(sort_cgmap_wrapper,l_multi_args)
    '''

    for str_cgmap in l_cgmap_files:
       sort_cgmap(str_cgmap,STR_TMP_DIR+str_out.split('.')[0],str_context,n_min_cov,n_sorted)

    #perform union
    write_intersection(l_sorted_cgmap,str_out)

#------------------------------------------------------------------------------
# main
#------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    main(args.in_sam,args.min_cov,args.out_sam,args.filter_context,args.sorted,args.proc_sam)