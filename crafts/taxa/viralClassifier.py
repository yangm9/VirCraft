import sys
from ..general import utils
from ..identify.viralDetectors import VirDetectTools
from .viralCompare import EnviComp

class VirTaxa(VirDetectTools):
    '''
    Taxonomic assignment
    '''
    def __init__(self, fasta=None, outdir=None, threads=8):
        super().__init__(fasta, outdir)
        self.threads = str(threads)
    def genomadTaxa(self):
        cmd, gn_dir = self.genomad()
        cmd.extend([utils.selectENV('VC-General')])
        gn_taxa = f'{gn_dir}/{self.name}_votus_summary/{self.name}_votus_virus_summary.tsv'
        votu_taxa = f'{self.outdir}/{self.name}.votu.taxa.txt'
        cmd.extend(['genomad_to_taxa.pl', gn_taxa, votu_taxa, '\n'])
        return cmd
    def Classify(self, unrun=False):
        cmd = [self.envs]
        tmp_cmd = self.genomadTaxa()
        cmd.extend(tmp_cmd)
        Comp = EnviComp(
            fasta=orf_faa,
            outdir=self.outdir,
            threads=self.threads
        )
        cmd.extend(Comp.vcontact2(vcont_db=self.confDict['vContact2DB']))
        shell = f'{self.shell_dir}/{self.name}_classify.sh'
        utils.printSH(shell, cmd)
        results = 0 if unrun else utils.execute(shell)
        return results
