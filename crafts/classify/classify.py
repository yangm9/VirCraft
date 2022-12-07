import sys
from ..general import cmdExec,general
from ..identify.viridsop import VirScan

class VirTaxa(VirScan):
    '''
    Taxonomic assignment
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir,threads)
    def demovir(self):
        '''
        Classify the virus contig by Demovir software for a certain single group.
        '''
        demovir_db=self.confDict['DemovirDB']
        TrEMBL_viral_taxa=f'{demovir_db}/TrEMBL_viral_taxa.RDS'
        TrEMBL_viral_taxa_lnk=f'{self.outdir}/TrEMBL_viral_taxa.RDS'
        cmd=['ln -s',TrEMBL_viral_taxa,TrEMBL_viral_taxa_lnk,'\n']
        uniprot_trembl_viral=f'{demovir_db}/uniprot_trembl.viral.udb'
        uniprot_trembl_viral_lnk=f'{self.outdir}/uniprot_trembl.viral.udb'
        cmd.extend(['ln -s',uniprot_trembl_viral,uniprot_trembl_viral_lnk,'\n'])
        demovir=f'{sys.path[0]}/bin/demovir.*'
        cmd.extend(['cp',demovir,self.outdir,'\n'])
        votus=f'{self.outdir}/03.vOTUs/merged_virus_positive_nodup.fa'
        demovir=f'{self.outdir}/demovir.sh'
        cmd.extend([demovir,votus,self.threads,'\n'])
        return cmd
    def Classify(self):
        cmd=[self.envs]
        cmd.extend(self.demovir())
        shell=f'{self.outdir}/{self.name}_demovir.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results

