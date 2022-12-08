from ..general import cmdExec,general
from ..config.config import Seq

class GeneFunc(Seq):
    '''
    Gene function annotation.
    1) Predict genes from fasta file;
    2) Annotate the function of these genes.
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.threads=str(threads)
    def genePred(self):
        wkdir=f'{self.outdir}/0.prodigal'
        general.mkdir(wkdir)
        temp_orfs_faa=f'{wkdir}/temp.orfs.faa'
        orfs_faa=f'{wkdir}/{self.name}_votus.faa'
        orfs_ffn=f'{wkdir}/{self.name}_votus.ffn'
        temp_orfs_ffn=f'{wkdir}/temp.orfs.ffn'
        temp_txt=f'{wkdir}/temp.txt'
        cmd=['prodigal','-i',self.fasta,'-a',temp_orfs_faa,
            '-d',temp_orfs_ffn,'-m','-o',temp_txt,'-p meta -q\n',
            'cut -f 1 -d \" \"',temp_orfs_faa,'>',orfs_faa,'\n',
            'cut -f 1 -d \" \"',temp_orfs_ffn,'>',orfs_ffn,'\n',
            f'rm -f {wkdir}/temp.*\n']
        return cmd,orfs_faa
    def eggnogAnno(self,orfs_faa):
        wkdir=f'{self.outdir}/1.eggnog'
        general.mkdir(wkdir)
        anno_prefix=f'{wkdir}/all_votus'
        seed_orth=f'{anno_prefix}.emapper.seed_orthologs'
        eggout=f'{wkdir}/eggout'
        eggnog_db=self.confDict['EggNOGDB']
        cmd=['emapper.py -m diamond --no_annot --no_file_comments'
            '--cpu',self.threads,'-i',orfs_faa,
            '-o',anno_prefix,'--data_dir',eggnog_db,'\n',
            'emapper.py','--annotate_hits_table',seed_orth,
            '--no_file_comments','-o',eggout,'--cpu',self.threads,
            '--data_dir',eggnog_db,'--override\n']
        return cmd
    def keggAnno(self,orfs_faa):
        wkdir=f'{self.outdir}/2.kegg'
        general.mkdir(wkdir)
        ko_prof=f'{self.confDict["KofamscanDB"]}/profiles'
        ko_list=f'{self.confDict["KofamscanDB"]}/ko_list'
        exec_anno=f'{wkdir}/all_votus.exec_annotation.txt'
        exec_anno_detail=f'{exec_anno}.xls'
        cmd=['exec_annotation -f mapper','--cpu',self.threads,
            '-p',ko_prof,'-k',ko_list,'-o',exec_anno,orfs_faa,'\n',
            'exec_annotation -f detail --cpu 32','-p',ko_prof,
            '-k',ko_list,'-o',exec_anno_detail,orfs_faa,'\n']
        return cmd
    def FuncAnnot(self):
        cmd=[self.envs]
        tmp_cmd,orfs_faa=self.genePred()
        cmd.extend(tmp_cmd)
        cmd.extend(self.eggnogAnno(orfs_faa))
        cmd.extend(self.keggAnno(orfs_faa))
        shell=f'{self.outdir}/{self.name}_gene_annot.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
