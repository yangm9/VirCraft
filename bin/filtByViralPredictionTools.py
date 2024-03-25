#!/usr/bin/env python3
import os,sys
import pandas as pd

BMVPTools=['deepvirfinder','virsorter','virsorter2','marvel','vibrant']
ToolsNum=len(BMVPTools)

def getWtPRes(ctg_l,wtpresults_d):
    WtPSet=set(ctg_l)
    for tool in BMVPTools:
        PosiCtgfile=f'{wtpresults_d}/identified_contigs_by_tools/{tool}.txt'
        cmd='cat '+PosiCtgfile+'|sort -u'
        tmp=set(os.popen(cmd).read().split('\n'))
        WtPSet=WtPSet.intersection(tmp)
    return WtPSet

def filtViruScaf(wtpresults_d,checkv_d):
    '''
    Scaffolds without a determined completeness and viral specific genes predicted by CheckV must contain viral signature using benchmarked viral prediction tools (DeepVirFinder, VirSorter, VirSorter2, MARVEL and VIBRANT) with conservative cutoff, published in standard operating procedure.
    '''
    checkv_qual_file=f'{checkv_d}/quality_summary.tsv'
    df=pd.read_csv(checkv_qual_file,sep='\t')
    CtgList=df['contig_id'].tolist()
    PosiCtgList=df[~((df['checkv_quality']=='Not-determined') & (df['viral_genes']==0))]['contig_id'].tolist()
    AmbiCtgList=df[(df['checkv_quality']=='Not-determined') & (df['viral_genes']==0)]['contig_id'].tolist()
    WtPSet=getWtPRes(CtgList,wtpresults_d)
    PosiCtgList.extend(list(WtPSet.intersection(set(AmbiCtgList))))
    for ctg in PosiCtgList: print(ctg)
    return 0

if __name__=='__main__':
    if len(sys.argv)==3:
        filtViruScaf(sys.argv[1],sys.argv[2])
    else:
        print(f'python {sys.argv[0]} <wtpresults_dir> <checkv_dir>')
