#!/usr/bin/env python3
import os
import sys
import re
import pandas as pd
import linkTab

#author: yangm@idsse.ac.cn
#version: 0.0.1 2023-04-12 14:51
#version: 0.0.2 2023-07-27 22:33

#Add the 1st columns name for a dramv annotation file
def add_1st_column_name(anno_f,column_1st_name,anno_c1_f):
    with open(anno_f) as f:
        lines=f.readlines()
    if len(lines)>0:
        lines[0]=column_1st_name+lines[0]
    else:
        print(f'Error: {anno_f} is an empty file')
    with open(anno_c1_f,'w') as f:
        f.writelines(lines)
    return 0

#Get the dictionary of consensus protein id for DRAM-v and VIBRANT
#e.g., geneid:contig_id:start±end
def get_common_gene_id_from_gff(gff,tool):
    geneidDict={}
    GFF=open(gff)
    for line in GFF:
        if line.startswith('#'):continue
        line=line.strip()
        gene_id=comm_gene_id=''
        items=line.split('\t')
        contig_id=items[0]
        if tool=='dram':
            gene_id=items[8].split(';')[0].replace('ID=','')
        elif tool=='vibrant':
            gene_id=items[8].split(';')[0].split('_')[1]
            gene_id=f'{contig_id}_{gene_id}'
        else:
            print('Error: tool must be dram or vibrant')
        comm_gene_id=f'{contig_id}:{items[3]}{items[6]}{items[4]}'
        geneidDict[gene_id]=comm_gene_id
    GFF.close()
    return geneidDict

#Extract the consensus protein id from faa file
def get_common_gene_id_for_vibrant(faa):
    geneidDict={}
    FAA=open(faa)
    for line in FAA:
        if not line.startswith('>'):continue
        line=line.strip()
        line=line.lstrip('>')
        items=line.split('\t')
        gene_id=items[0]
        contig_id=gene_id.rsplit('_',1)[0]
        start,end=items[1].split('..')
        start=start.replace('(','')
        end=end.replace(')','')
        if items[2]=='1':
            strand='+'
        else:
            strand='-'
        comm_gene_id=f'{contig_id}:{start}{strand}{end}'
        geneidDict[gene_id]=comm_gene_id
    FAA.close()
    return geneidDict

#Output the private protein ID and its corresponding consensus protein ID
def list_common_gene_id(dm_gff,vb_faa,wkdir):
    dm_geneidDict=get_common_gene_id_from_gff(dm_gff,'dram')
    #vb_geneidDict=get_common_gene_id_from_gff(vb_gff,'vibrant')
    vb_geneidDict=get_common_gene_id_for_vibrant(vb_faa)
    dm_gene_id_list=f'{wkdir}/dram_id.tsv'
    DMGI=open(dm_gene_id_list,'w')
    DMGI.write('protein_id\tcommon_protein_id\n')
    for gene_id in dm_geneidDict.keys():
        line=f'{gene_id}\t{dm_geneidDict[gene_id]}\n'
        DMGI.write(line)
    DMGI.close()
    vb_gene_id_list=f'{wkdir}/vibrant_id.tsv'
    VBGI=open(vb_gene_id_list,'w')
    VBGI.write('protein\tcommon_protein_id\n')
    for gene_id in vb_geneidDict.keys():
        pattern=r'_fragment_\d+'
        common_prot_id=re.sub(pattern,'',vb_geneidDict[gene_id])
        line=f'{gene_id}\t{common_prot_id}\n'
        VBGI.write(line)
    VBGI.close()
    return dm_gene_id_list,vb_gene_id_list

def combine_amgs_of_dramv_and_vibrant(dmdir,vbdir,wkdir='.'):
    vbdir.strip('/')
    seq_name=os.path.basename(vbdir).split('_',1)[1]
    vb_faa=f'{vbdir}/VIBRANT_phages_{seq_name}/{seq_name}.phages_combined.faa'
    #vb_gff=f'{vbdir}/{seq_name}.prodigal.gff'
    dm_gff=f'{dmdir}/genes.gff'
    #Get the common gene ID list of dramv and vibrant base on respective gff.
    dm_gene_id_list,vb_gene_id_list=list_common_gene_id(dm_gff,vb_faa,wkdir)
    dm_annot=f'{dmdir}/annotations.tsv'
    dm_annot_c1=f'{wkdir}/annotations_c1.tsv'
    #add the 1st column name for dramv
    add_1st_column_name(dm_annot,'protein_id',dm_annot_c1)
    vb_annot=f'{vbdir}/VIBRANT_results_{seq_name}/VIBRANT_annotations_{seq_name}.tsv'
    #Add the common gene id for dramv and vibrant annotation file, respectively
    dm_id_annot=f'{wkdir}/dram_annotations.tsv'
    linkTab.merge(dm_annot_c1,dm_gene_id_list,'left','protein_id',dm_id_annot)
    vb_id_annot=f'{wkdir}/vibrant_annotations.tsv'
    linkTab.merge(vb_annot,vb_gene_id_list,'left','protein',vb_id_annot)
    #Merge the annotation file of dramv and vibrant by common gene id
    merged_annot=f'{wkdir}/merged_annotations.tsv'
    try:
        linkTab.merge(dm_id_annot,vb_id_annot,'outer','common_protein_id',merged_annot)
    except Exception as e:
        print('ERROR: No Match!!!')
        exit(1)
    return merged_annot

#Filter AMGs according to the rules provided by Zhou et al.(10.1038/s42003-022-04027-y)
def generate_amgs_table(dm_vb_merged_annot):
    rename_columns_dict={
        'common_protein_id':'gene_id','protein_id':'dram_protein_id',
        'scaffold_x':'dram_scaffold','gene_position':'gene_position',
        'start_position':'start','end_position':'end','strandedness':'strand',
        'rank':'dram_rank','protein':'vibrant_protein_id',
        'scaffold_y':'vibrant_scaffold','KO name':'KO_name',
        'KO evalue':'KO_evalue','KO score':'KO_score',
        'KO v-score':'KO_v-score','Pfam name':'Pfam_name',
        'Pfam evalue':'Pfam_evalue','Pfam score':'Pfam_score',
        'Pfam v-score':'Pfam_v-score','VOG name':'VOG_name',
        'VOG evalue':'VOG_evalue','VOG score':'VOG_score',
        'VOG v-score':'VOG_v-score'
    }
    amg_criteria="(not (auxiliary_score>3 or 'T' in amg_flags or 'B' in amg_flags)) and AMG=='AMG'"
    df=pd.read_csv(dm_vb_merged_annot,sep='\t',low_memory=False)
    #Put the common_protein_id column to 1st
    cols=df.columns.tolist()
    cols.remove('common_protein_id')
    cols.insert(0,'common_protein_id')
    df=df.reindex(columns=cols)
    df=df.drop('fasta', axis=1)
    df=df.rename(columns=rename_columns_dict)
    wkdir=os.path.dirname(dm_vb_merged_annot)
    all_merged_annot=f'{wkdir}/all_merged_annotation.tsv'
    df.to_csv(all_merged_annot,sep='\t',index=False)
    df=df.query(amg_criteria)
    filt_merged_annot=f'{wkdir}/filted_merged_amg.tsv'
    df.to_csv(filt_merged_annot,sep='\t',index=False)
    return 0

if __name__=='__main__':
    dm_vb_merged_annot=combine_amgs_of_dramv_and_vibrant(sys.argv[1],sys.argv[2],sys.argv[3])
    generate_amgs_table(dm_vb_merged_annot)
    os.remove(dm_vb_merged_annot)
