import argparse
from rtree import index
import re

''' demo index
    idx = index.Index()
    start, chrid, end, chrid = (0.0, 0.0, 1.0, 0.0)
    idx.insert(0, (start, chrid, end, chrid))
    lval = list(idx.intersection((1.1, 0.0, 1.1, 0.0)))
'''

def test_bbox():
    idx = index.Index()
    start, chrid, end, chrid = (0.0, 0.0, 1.0, 0.0)
    idx.insert(0, (start, chrid, end, chrid))
    lhit = idx.intersection((1.1, 0.0, 1.1, 0.0),objects=True)
    for i in lhit:
        nstart,nchr,nend,nchr2 = i.bbox
        print '%d:%d-%d'%(nchr,nstart,nend)

    '''
    if len(lval) > 0:
        print '%d'%lval[0]
    '''

def chrlbl_to_num(str_chr,l_chr_lbl):
    n_chr_id = 0
    if str_chr in l_chr_lbl:
        n_chr_id = l_chr_lbl.index(str_chr)
    else:
        l_chr_lbl.append(str_chr)
        n_chr_id = len(l_chr_lbl)
    return n_chr_id

def num_to_chrlbl(nnum,l_chr_lbl):
    if nnum >= len(l_chr_lbl):
        print 'chr label index out of range %d'%nnum
    return l_chr_lbl[nnum]

def load_mappable(str_mappable_file,l_chr_lbl):
    frag_idx = index.Index()
    for str_line in open(str_mappable_file,'r'):
        lcols = str_line.rstrip().split('\t')
        n_chr = chrlbl_to_num(lcols[0],l_chr_lbl)
        n_frag_id = int(lcols[1])
        n_start = int(lcols[2])+1
        n_end = int(lcols[3])+1
        frag_idx.insert(n_frag_id,(n_start,n_chr,n_end,n_chr))
    return frag_idx

def load_pos_fcgmap(str_comm_cgmap,l_chr_lbl):
    l_pos = []
    fin = open(str_comm_cgmap,'r')
    b_is_header = True
    for str_line in fin:
        if b_is_header:
            b_is_header = False
            continue
        lcols = str_line.rstrip().split('\t')
        l_chr_lbls = lcols[0].split('_')
        str_chr = '_'.join(l_chr_lbls[:-1])
        str_pos = l_chr_lbls[-1]
        n_chr = chrlbl_to_num(str_chr,l_chr_lbl)
        n_pos = int(str_pos)
        l_pos.append((n_pos,n_chr,n_pos,n_chr))
    fin.close()
    return l_pos

def load_comm_pos(str_pos,l_chr_lbl):
    l_pos = []
    for str_line in open(str_pos):
        lcols = str_line.rstrip().split(' ')
        n_chr = chrlbl_to_num(lcols[0],l_chr_lbl)
        n_pos = int(lcols[1])
        l_pos.append((n_pos,n_chr,n_pos,n_chr))
    return l_pos

def get_pos_to_frag(frag_idx,lpos,l_chr_lbl):
    d_pos_to_frag = {}
    for i in range(len(lpos)):
        l_cur_pos = lpos[i]

        str_cur_pos = '%s_%d'%(num_to_chrlbl(l_cur_pos[1],l_chr_lbl),l_cur_pos[0])
        l_hits = frag_idx.intersection(l_cur_pos,objects=True)
        for i in l_hits:
            nstart,nchr,nend,nchr2 = i.bbox
            str_frag = '%s:%d-%d'%(num_to_chrlbl(l_cur_pos[1],l_chr_lbl),int(nstart),
                                   int(nend))
            d_pos_to_frag[str_cur_pos] = str_frag
            break
    return d_pos_to_frag

def proc_fcgmap(d_pos_to_frag,str_fcgmap,str_out):
    b_header = True
    fout = open(str_out,'w')
    for str_line in open(str_fcgmap,'r'):
        lcols = str_line.rstrip().split('\t')
        if b_header:
            fout.write('\t'.join(lcols[0:1])+'\tFRAG\t'+'\t'.join(lcols[1:])+'\n')
            b_header = False
        else:
            str_frag = d_pos_to_frag[lcols[0]]
            fout.write('\t'.join(lcols[0:1])+'\t'+str_frag+'\t'+'\t'.join(lcols[1:])+'\n')
    fout.close()


#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pos', type=str, help="comm position file",default='')
    parser.add_argument('-m', '--mappable', type=str, help="mappable regions file")
    parser.add_argument('-f', '--fcgmap', type=str, help="sample methylation matrix")
    args = parser.parse_args()

    #list of chr labels
    l_chr_lbls = []

    #load fragments into rtree
    frag_idx = load_mappable(args.mappable,l_chr_lbls)

    lpos = []
    if args.pos == '':
        lpos = load_pos_fcgmap(args.fcgmap,l_chr_lbls)
    else:
        lpos = load_comm_pos(args.pos)

    d_pos_to_frag = get_pos_to_frag(frag_idx,lpos,l_chr_lbls)

    proc_fcgmap(d_pos_to_frag,args.fcgmap,re.sub('.txt$','_frag.txt',args.fcgmap))

    print'end'

