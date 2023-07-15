import os
from datetime import datetime
from ..general import utils

class DB:
    '''
    Deploy the databases of VirCraft
    '''
    def __init__(self,outdir='',threads=8):
        self.outdir=outdir
        self.threads=str(threads)
    def dl_refseq_viral_prot(self):
        wkdir=f'{self.outdir}/tax_db_dir'
        utils.mkdir(wkdir)
        vir_prot_url='https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz'
        vir_prot=f'{wkdir}/viral.1.protein.faa.gz'
        cmd=['wget -c',vir_prot_url,'-O',vir_prot,'--no-check-certificate\n',
            'gzip -d',vir_prot,'\n']
        return cmd
    def dl_viridsop_db(self):
        wkdir=f'{self.outdir}/vs2_db'
        utils.mkdir(wkdir)
        cmd=['conda run -n viral-id-sop virsorter setup',
             '-d',wkdir,'-j',threads,'\n']
        return cmd
    def dl_checkv_db(self):
        wkdir=f'{self.outdir}'
        utils.mkdir(wkdir)
        cmd=['conda run -n viral-id-sop checkv download_database',wkdir,'\n']
        return cmd
    
    def dl_dramv_db(self):
        wkdir=f'{self.outdir}/dramv_db'
        utils.mkdir(wkdir)
        cmd=['conda run -n viral-id-sop DRAM-setup.py prepare_databases',
             '--skip_uniref','--output_dir',wkdir,'\n']
        return cmd
    def dl_vb_db(self):
        wkdir=f'{self.outdir}/vb_db'
        utils.mkdir(wkdir)
        cmd=['conda run -n vibrant download-db.sh',wkdir,'\n']
        return cmd
#    def 
#    def dl_all(self):
        


