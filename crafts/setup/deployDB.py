import os
import sys
import shutil
from datetime import datetime
from ..general import utils
from . import URL

class DB:
    '''
    Deploy the databases of VirCraft
    '''
    def __init__(self,outdir='',threads=8):
        self.outdir=os.path.abspath(outdir)
        self.threads=str(threads)
    def dl_virsorter2_db(self):
        wkdir=f'{self.outdir}/VC-VirSorter2DB'
        utils.mkdir(wkdir)
        cmd=['conda run -n VC-VirSorter2 virsorter setup',
             '-d',wkdir,'-j',self.threads,'\n']
        return cmd
    def dl_vibrant_db(self):
        wkdir=f'{self.outdir}/VC-VIBRANTDB'
        utils.mkdir(wkdir)
        cmd=['conda run -n VC-VIBRANT download-db.sh',wkdir,'\n']
        return cmd
#    def dl_deepvirfinder_db(self):
#        wkdir=f'{self.outdir}/VC-DeepVirFinder'
#        utils.mkdir(wkdir)
#        which_dvf_cmd='conda run -n VC-DeepVirFinder which dvf.py'
#        models_dir=os.path.dirname(os.popen(which_dvf_cmd).read().strip())
#        models_dir=models_dir.replace('VC-DeepVirFinder/bin','VC-DeepVirFinder/share/deepvirfinder/models')
#        cmd=['cp -r',models_dir,wkdir,'\n']
#        return cmd
    def dl_CAT_db(self):
        wkdir=f'{self.outdir}/VC-CATDB'
        db_basename=os.path.basename(URL.CAT_DB_URL)
        db_file=f'{wkdir}/{db_file_basename}'
        cmd=['wget -c',URL.CAT_DB_URL,'-O',db_file,'\n',
            '&& tar xzf',db_file,'-C',wkdir,'\n']
        return cmd
    def dl_checkv_db(self):
        wkdir=f'{self.outdir}/VC-CheckVDB'
        db_dir=f'{wkdir}/checkv-db-v*'
        db_files=db_dir+'/*'
        utils.mkdir(wkdir)
        cmd=['conda run -n VC-CheckV checkv download_database',wkdir,'\n'
            'mv',db_files,wkdir,'&& rmdir',db_dir,'\n']
        return cmd
    def dl_refseq_viral_prot(self):
        wkdir=f'{self.outdir}/VC-ViralRefSeqDB'
        utils.mkdir(wkdir)
        taxdump_tgz=f'{wkdir}/taxdump.tar.gz'
        vir_prot_prefix=f'{wkdir}/viral.1.protein'
        vir_prot=vir_prot_prefix+'.faa'
        vir_prot_gz=vir_prot+'.gz'
        vir_version=f'{wkdir}/RELEASE_NUMBER'
        vir_pid_sp=f'{wkdir}/NCBI_viral_pid_sp.txt'
        vir_name_taxaid=vir_pid_sp.replace('_pid_sp','_name_taxid')
        vir_pid_taxa=vir_pid_sp.replace('_pid_sp','_taxnomomy')
        vir_sp_taxa=vir_pid_sp.replace('_pid_sp','_sp_taxa')
        vir_pid_sp_h=vir_pid_sp.replace('.txt','.h.txt')
        vir_full_taxa=vir_pid_sp.replace('_pid_sp','_full_taxnomomy')
        cmd=[utils.selectENV('VC-General')]
        cmd.extend(
            ['wget','-c',URL.NCBI_VIR_PROT_URL,'-O',vir_prot_gz,
            '--no-check-certificate\n','gzip -d',vir_prot_gz,'\n',
            'wget','-c',URL.NCBI_RELEASE_NUMBER_URL,'-O',vir_version,
            '--no-check-certificate\n','makeblastdb -in',vir_prot,
            '-parse_seqids -hash_index','-out',vir_prot_prefix,'-dbtype prot\n',
            'wget','-c',URL.NCBI_TAXDUMP_URL,'-O',taxdump_tgz,'--no-check-certificate\n',
            '''if [ ! -d ~/.taxonkit ]; then
    mkdir ~/.taxonkit
fi
if [ -z "$(ls -A ~/.taxonkit)" ]; then
    tar xzf''',taxdump_tgz,'''-C ~/.taxonkit
fi
'''
            'extract_pid_sp_from_faa.py',vir_prot,'>',vir_pid_sp,'\n',
            'cut -f 2',vir_pid_sp,"|uniq|taxonkit name2taxid|sed '1ispecies\\ttaxid'>",vir_name_taxaid,'\n',
            'cut -f 2',vir_pid_sp,"|uniq|taxonkit name2taxid|cut -f 2|taxonkit lineage|taxonkit reformat -r 'Unassigned'|cut -f 1,3|sed '1itaxid\\ttaxonomy'>",vir_pid_taxa,'\n',
            "csvtk join -t -f 'taxid;taxid'",vir_name_taxaid,vir_pid_taxa,
            '|uniq>',vir_sp_taxa,'\n',
            'csvtk add-header -t -n NCBI_ID,species',vir_pid_sp,
            '>',vir_pid_sp_h,'\n',
            "csvtk join -t -f 'species;species'",vir_sp_taxa,vir_pid_sp_h,
            '>',vir_full_taxa,'\n',
            'rm -f',taxdump_tgz,vir_pid_sp,vir_name_taxaid,
            vir_pid_taxa,vir_sp_taxa,vir_pid_sp_h,'\n']
        )
        return cmd
    def dl_demovir_db(self):
        wkdir=f'{self.outdir}/VC-DemovirDB'
        utils.mkdir(wkdir)
        TrEMBL_viral_taxa_rds=f'{wkdir}/TrEMBL_viral_taxa.RDS'
        nr95_bz2=f'{wkdir}/nr.95.fasta.bz2'
        uniprot_udb=f'{wkdir}/uniprot_trembl.viral.udb'
        nr95_fa=nr95_bz2.replace('nr.95.fasta.bz2','nr.95.fasta')
        cmd=[utils.selectENV('VC-General')]
        cmd.extend(
            ['wget','-c',URL.Demovir_URL,'-O',TrEMBL_viral_taxa_rds,'\n',
            'wget','-c',URL.DemovirDB_URL,'-O',nr95_bz2,
            '&& cd',wkdir,'&& bzip2 -d',nr95_bz2,'\n',
            'usearch','-makeudb_ublast',nr95_fa,'-output',uniprot_udb,
            '> usearch_database.log\n']
        )
        return cmd
    def dl_eggnog_db(self):
        wkdir=f'{self.outdir}/VC-eggNOGDB'
        utils.mkdir(wkdir)
        eggnog_db_gz=f'{wkdir}/eggnog.db.gz'
        eggnog_dmnd_gz=f'{wkdir}/eggnog_proteins.dmnd.gz'
        cmd=['wget -c',URL.EGGNOG_DB_URL,'-O',eggnog_db_gz,
        '&& gunzip',eggnog_db_gz,'\n',
        'wget -c',URL.EGGNOG_DMND_URL,'-O',eggnog_dmnd_gz,
        '&& gunzip',eggnog_dmnd_gz,'\n']
        return cmd
    def dl_kofamscan_db(self):
        wkdir=f'{self.outdir}/VC-KofamScanDB'
        utils.mkdir(wkdir)
        ko_list_gz=f'{wkdir}/ko_list.gz'
        profiles_tgz=f'{wkdir}/profiles.tar.gz'
        cmd=['wget -c',URL.KO_LIST_URL,'-O',ko_list_gz,'\n',
            'gunzip',ko_list_gz,'\n',
            'wget -c',URL.KO_PROFILES_URL,'-O',profiles_tgz,'\n',
            'tar xzf',profiles_tgz,'-C',wkdir,'\n']
        return cmd
    def dl_dramv_db(self):
        wkdir=f'{self.outdir}/VC-DRAMvDB'
        utils.mkdir(wkdir)
        cmd=['conda run -n VC-DRAMv DRAM-setup.py prepare_databases',
             '--skip_uniref','--output_dir',wkdir,'\n']
        return cmd
    def dl_gtdbtk_db(self):
        wkdir=f'{self.outdir}/VC-GTDBTkDB'
        utils.mkdir(wkdir)
        gtdbtk_db_gz=os.path.basename(URL.GTDBTK_DB_URL)
        virsion=gtdbtk_db_gz.replace('_data.tar.gz','')
        virsion=virsion.replace('gtdbtk_r','')
        db_virsion_dir=wkdir+'/release'+virsion
        dbdir=wkdir+'/GTDBTkDB'
        gtdbtk_db_gz=os.path.join(wkdir,gtdbtk_db_gz)
        cmd=['wget -c',URL.GTDBTK_DB_URL,'-O',gtdbtk_db_gz,
        '--no-check-certificate\n',
        'tar xzf',gtdbtk_db_gz,'-C',wkdir,'\n',
        'mv',db_virsion_dir,dbdir,'\n',
        'rm -f',gtdbtk_db_gz,'\n']
        gtdbtk_data_path=f'GTDBTK_DATA_PATH={dbdir}'
        cmd.extend(
            ['conda env config vars set',gtdbtk_data_path,'-n VC-GTDBTk\n']
        )
        return cmd
    def Deploy(self,unrun=False,clear=False):
        cmd=self.dl_virsorter2_db()
        shell=f'{self.outdir}/virsorter2_db_deploy.sh'
        utils.printSH(shell,cmd)
        cmd=self.dl_vibrant_db()
        shell=f'{self.outdir}/vibrant_db_deploy.sh'
        utils.printSH(shell,cmd)
 #       cmd=self.dl_deepvirfinder_db()
 #       shell=f'{self.outdir}/deepvirfinder_db_deploy.sh'
 #       utils.printSH(shell,cmd) 
        cmd=self.dl_checkv_db()
        shell=f'{self.outdir}/CAT_db_deploy.sh'
        utils.printSH(shell,cmd)
        cmd=self.dl_checkv_db()
        shell=f'{self.outdir}/checkv_db_deploy.sh'
        utils.printSH(shell,cmd)
        cmd=self.dl_refseq_viral_prot()
        shell=f'{self.outdir}/ncbirefseq_db_deploy.sh'
        utils.printSH(shell,cmd)
        cmd=self.dl_demovir_db()
        shell=f'{self.outdir}/demovir_db_deploy.sh'
        utils.printSH(shell,cmd)
        cmd=self.dl_eggnog_db()
        shell=f'{self.outdir}/eggnog_db_deploy.sh'
        utils.printSH(shell,cmd)
        cmd=self.dl_kofamscan_db()
        shell=f'{self.outdir}/kegg_db_deploy.sh'
        utils.printSH(shell,cmd)
        cmd=self.dl_dramv_db()
        shell=f'{self.outdir}/dramv_db_deploy.sh'
        utils.printSH(shell,cmd)
        cmd=self.dl_gtdbtk_db()
        shell=f'{self.outdir}/gtdbtk_db_deploy.sh'
        utils.printSH(shell,cmd)
        cmd=[utils.selectENV('VC-VIBRANT')]
        sed_cmd=f"sed 's/\/data_backup\/database/{self.outdir}/'"
        config_template=f'{sys.path[0]}/crafts/config/config.tpl'
        config_file=f'{sys.path[0]}/config'
        cmd.extend(
            ['multithreads.pl',self.outdir,'db_deploy.sh 2\n',
            sed_cmd,config_template,'>',config_file,'\n']
        )
        shell=f'{self.outdir}/all_database_deploy.sh'
        utils.printSH(shell,cmd)
        results=''
        if not unrun: results=utils.execute(shell)
        return results
