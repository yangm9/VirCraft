import os
import sys
import pandas as pd
import linkTab

def getCommGeneID(gff,tool):
    geneidDict={}
    GFF=open(gff)
    for line in GFF:
        if line.startswith('#'):
            continue
        line=line.strip()
        contig_id=gene_id=comm_gene_id=''
        items=line.split('\t')
        if tool=='dram':
            contig_id=items[0].split('__')[0]
            gene_id=items[8].split(';')[0].replace('ID=','')
        elif tool=='vibrant':
            contig_id=items[0]
            gene_id=items[8].split(';')[0].split('_')[1]
            gene_id=f'{contig_id}_{gene_id}'
        else:
            print('Error: tool must be dram or vibrant')
        comm_gene_id=f'{contig_id}_{items[3]}_{items[4]}'
        geneidDict[gene_id]=comm_gene_id
    GFF.close()
    return geneidDict

def list_common_gene_id(dm_gff,vb_gff,wkdir):
    dm_geneidDict=getCommGeneID(dm_gff,'dram')
    vb_geneidDict=getCommGeneID(vb_gff,'vibrant')
    dm_gene_id_list=f'{wkdir}/dram_id.li'
    DMGI=open(dm_gene_id_list,'w')
    DMGI.write('protein_id\tcommon_protein_id\n')
    for gene_id in dm_geneidDict.keys():
        line=f'{gene_id}\t{dm_geneidDict[gene_id]}\n'
        DMGI.write(line)
    DMGI.close()
    vb_gene_id_list=f'{wkdir}/vibrant_id.li'
    VBGI=open(vb_gene_id_list,'w')
    VBGI.write('protein\tcommon_protein_id\n')
    for gene_id in vb_geneidDict.keys():
        line=f'{gene_id}\t{vb_geneidDict[gene_id]}\n'
        VBGI.write(line)
    VBGI.close()
    return dm_gene_id_list,vb_gene_id_list

def add_1st_column_name(anno_f,column_1st_name,anno_c1_f):
    "Add the 1st columns name for dramv annotation file"
    with open(anno_f) as f:
        lines=f.readlines()
    if len(lines)>0:
        lines[0]=column_1st_name+lines[0]
    else:
        print(f'Error: {anno_f} is an empty file')
    with open(anno_c1_f,'w') as f:
        f.writelines(lines)
    return 0

def combine_amgs_of_dramv_and_vibrant(dmdir,vbdir,wkdir='.'):
    vbdir.strip('/')
    seq_name=os.path.basename(vbdir).split('_',1)[1]
    vb_gff=f'{vbdir}/{seq_name}.prodigal.gff'
    dm_gff=f'{dmdir}/genes.gff'
    #Get the common gene ID list of dramv and vibrant base on respective gff.
    dm_gene_id_list,vb_gene_id_list=list_common_gene_id(dm_gff,vb_gff,wkdir)
    dm_annot=f'{dmdir}/annotations.tsv'
    dm_annot_c1=f'{wkdir}/annotations_c1.xls'
    #add the 1st column name for dramv
    add_1st_column_name(dm_annot,'protein_id',dm_annot_c1)
    vb_annot=f'{vbdir}/VIBRANT_results_{seq_name}/VIBRANT_annotations_{seq_name}.tsv'
    #Add the common gene id for dramv and vibrant annotation file, respectively
    dm_id_annot=f'{wkdir}/dram_annotations.xls'
    linkTab.merge(dm_annot_c1,dm_gene_id_list,'left','protein_id',dm_id_annot)
    vb_id_annot=f'{wkdir}/vibrant_annotations.xls'
    linkTab.merge(vb_annot,vb_gene_id_list,'left','protein',vb_id_annot)
    #Merge thw annotation file of dramv and vibrant by common gene id
    merged_annot=f'{wkdir}/merged_annotations.xls'
    linkTab.merge(dm_id_annot,vb_id_annot,'outer','common_protein_id',merged_annot)
    return merged_annot

#def amgs_filter(dm_vb_merged_annot):

if __name__=='__main__':
    combine_amgs_of_dramv_and_vibrant(sys.argv[1],sys.argv[2],sys.argv[3])
