from datetime import datetime
from ..general import utils
from . import URL

class DB:
    '''
    Deploy the databases of VirCraft
    '''
    def __init__(self,outdir='',threads=8):
        self.outdir=outdir
        self.threads=str(threads)
    def dl_virsorter2_db(self):
        wkdir=f'{self.outdir}/VC-VirSorter2DB'
        utils.mkdir(wkdir)
        cmd=['conda run -n viral-id-sop virsorter setup',
             '-d',wkdir,'-j',self.threads,'\n']
        return cmd
    def dl_vibrant_db(self):
        wkdir=f'{self.outdir}/VC-VIBRANTDB'
        utils.mkdir(wkdir)
        cmd=['conda run -n vibrant download-db.sh',wkdir,'\n']
        return cmd
    def dl_checkv_db(self):
        wkdir=f'{self.outdir}/VC-CheckVDB'
        utils.mkdir(wkdir)
        cmd=['conda run -n viral-id-sop checkv download_database',wkdir,'\n']
        return cmd
    def dl_refseq_viral_prot(self):
        wkdir=f'{self.outdir}/VC-ViralRefSeqDB'
        utils.mkdir(wkdir)
        taxdump_tgz=f'{wkdir}/taxdump.tar.gz'
        vir_prot_prefix=f'{wkdir}/viral.1.protein'
        vir_prot=vir_prot_prefix+'.faa'
        vir_prot_gz=vir_prot+'.gz'
        vir_pid_sp=f'{wkdir}/NCBI_viral_pid_sp.txt'
        vir_name_taxaid=vir_pid_sp.replace('_pid_sp','_name_taxid')
        vir_pid_taxa=vir_pid_sp.replace('_pid_sp','_taxnomomy')
        vir_sp_taxa=vir_pid_sp.replace('_pid_sp','_sp_taxa')
        vir_pid_sp_h=vir_pid_sp.replace('.txt','.h.txt')
        vir_full_taxa=vir_pid_sp.replace('_pid_sp','_full_taxnomomy')
        cmd=['wget','-c',URL.NCBI_VIR_PROT_URL,'-O',vir_prot_gz,
            '--no-check-certificate\n','gzip -d',vir_prot_gz,'\n',
            'wget','-c',URL.NCBI_RELEASE_NUMBER_URL,'-O',vir_prot,
            '--no-check-certificate\n','makeblastdb -in',vir_prot_gunzip,
            '-parse_seqids -hash_index','-out',vir_prot_prefix,'-dbtype prot\n',
            'wget','-c',URL.NCBI_TAXDUMP_URL,'-O',taxdump_tgz,'--no-check-certificate\n',
            'mkdir ~/.taxonkit && tar xzf',taxdump_tgz,'-c ~/.taxonkit\n',
            'extract_name.py',vir_prot,'>',vir_pid_sp,'\n',
            'cut -f 2',vir_pid_sp,"|uniq|taxonkit name2taxid|sed '1ispecies\ttaxid'>",vir_name_taxaid,'\n',
            'cut -f 2',vir_pid_sp,"|uniq|taxonkit name2taxid|cut -f 2|taxonkit lineage|taxonkit reformat -r 'Unassigned'|cut -f 1,3|sed '1itaxid\ttaxonomy'>",vir_pid_taxa,'\n',
            "csvtk join -t -f 'taxid;taxid'",vir_name_taxaid,vir_pid_taxa,
            '|uniq>',vir_sp_taxa,'\n',
            'csvtk add-header -t -n NCBI_ID,species',vir_pid_sp,
            '>',vir_pid_sp_h,'\n',
            "csvtk join -t -f 'species;species'",vir_sp_taxa,vir_pid_sp_h,
            '>',vir_full_taxa,'\n',
            'rm -f',vir_pid_sp,vir_name_taxaid,vir_pid_taxa,vir_sp_taxa,vir_pid_sp_h,'\n']
        return cmd
    def dl_kofamscan_db(self):
        wkdir=f'{self.outdir}/VC-KofamScanDB'
        utils.mkdir(wkdir)
        ko_list_gz=f'{wkdir}/ko_list.gz'
        profiles_tgz=f'{wkdir}/profiles.tar.gz'
        cmd=['wget -c',URL.KO_LIST_URL,'-O',ko_list_gz,'\n',
            'gunzip',ko_list_gz,'\n',
            'wget -c',URL.KO_PROFILES_URL,'-O',profiles_gz,'\n',
            'tar xzf',profiles_tgz,'-C',wkdir,'\n']
        return cmd
    def dl_dramv_db(self):
        wkdir=f'{self.outdir}/VC-DRAMvDB'
        utils.mkdir(wkdir)
        cmd=['conda run -n viral-id-sop DRAM-setup.py prepare_databases',
             '--skip_uniref','--output_dir',wkdir,'\n']
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
    def Deploy(self,unrun=False,clear=False):
        cmd=self.dl_virsorter2_db()
        shell=f'{self.outdir}/virsorter2_db_deploy.sh'
        utils.printSH(shell,cmd)
        cmd=self.dl_vibrant_db()
        shell=f'{self.outdir}/vibrant_db_deploy.sh'
        utils.printSH(shell,cmd)
        cmd=self.dl_checkv_db()
        shell=f'{self.outdir}/checkv_db_deploy.sh'
        utils.printSH(shell,cmd)
        cmd=self.dl_refseq_viral_prot()
        shell=f'{self.outdir}/ncbirefseq_db_deploy.sh'
        utils.printSH(shell,cmd)
        cmd=self.dl_dramv_db()
        shell=f'{self.outdir}/dramv_db_deploy.sh'
        utils.printSH(shell,cmd)
        cmd.extend(['multithreads.pl',self.outdir,'db_deploy.sh 2\n'])
        results=''
        if not unrun: results=utils.execute(cmd)
        return results