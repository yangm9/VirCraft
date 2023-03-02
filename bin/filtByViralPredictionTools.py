#!/usr/bin/env python3
import os,sys
import pandas as pd

def filtViruScaf(wtpresults_d,checkv_d):
    "Scaffolds without a determined completeness and viral specific genes predicted by CheckV must contain viral signature using benchmarked viral prediction tools (DeepVirFinder, VirSorter, VirSorter2, MARVEL and VIBRANT) with conservative cutoff, published in standard operating procedure"

    BMVPTools=['deepvirfinder','virsorter','virsorter2','marvel','vibrant']
    posi_scaf_files=wtpresults_d+'/identified_contigs_by_tools/*.txt'
    cmd='cat '+posi_scaf_files+'|sort -u'
    PosiToolScafList=os.popen(cmd).read().split('\n')
    complete_genomes_tsv=checkv_d+'/complete_genomes.tsv'
    df=pd.read_csv(complete_genomes_tsv,sep='\t')
    CompleteGenomeList=df['contig_id'].tolist()
    PosiScafList=list(set(PosiToolScafList+CompleteGenomeList))
    quality_summary_tsv=checkv_d+'/quality_summary.tsv'
    QST=open(quality_summary_tsv)
    head=QST.readline()
    for line in QST:
        items=line.strip('\n').split('\t')
        seq_name=items[0]
        if seq_name in PosiScafList or int(items[5])>0:
            print(seq_name)
    QST.close()
    return 0

if __name__=='__main__':
    if len(sys.argv)==3:
        filtViruScaf(sys.argv[1],sys.argv[2])
    else:
        print(f'python {sys.argv[0]} <wtpresults_dir> <checkv_dir>')
