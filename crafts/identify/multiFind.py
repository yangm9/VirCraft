import os
from ..general import utils
from ..identify.viridsop import VirScan

class MultiTools(VirScan):
    '''
    Generate command for DeepVirFinder, VIBRANT and Contig Annotation Tool (CAT).
    '''
    def __init__(self, fasta = '', outdir = '', threads = 8):
        super().__init__(fasta, outdir)
        self.threads = str(threads)
    def deepvirfinder(self, cutoff: str):
        cmd = [utils.selectENV('VC-DeepVirFinder')]
        wkdir = f'{self.outdir}/deepvirfinder'
        utils.mkdir(wkdir)
        cmd.extend(
            ['dvf.py', '-i', self.fasta, '-m', self.confDict['DeepVirFinderDB'],  '-o', wkdir, '-c', self.threads, '-l', cutoff, '\n']
        )
        return cmd, wkdir
    def vibrant(self, cutoff=1500):
        cutoff=str(cutoff)
        cmd = [utils.selectENV('VC-VIBRANT')]
        wkdir = f'{self.outdir}/VIBRANT_{self.name}'
        VIBRANT_DB_databases = self.confDict['VIBRANTDB'] + '/databases'
        VIBRANT_DB_files = self.confDict['VIBRANTDB'] + '/files'
        cmd.extend(
            ['VIBRANT_run.py', '-i', self.fasta, '-l', cutoff, '-t', self.threads, '-folder', self.outdir, '-d', VIBRANT_DB_databases, '-m', VIBRANT_DB_files, '\n']
        )
        return cmd, wkdir
    def genomad(self):
        cmd = [utils.selectENV('VC-GeNomad')]
        wkdir = f'{self.outdir}/genomad'
        GeNomad_DB = self.confDict['GeNomadDB']
        cmd.extend(
            ['genomad end-to-end', '--cleanup', '--threads', self.threads, self.fasta, wkdir, GeNomad_DB, '\n']
        )
        return cmd, wkdir
    def contig_annotation_tool(self, in_fa=None):
        cmd = [utils.selectENV('VC-General')]
        wkdir = f'{self.outdir}/CAT_{self.name}'
        utils.mkdir(wkdir)
        CAT_DB = self.confDict['CATDB'] + '/db'
        CAT_TAX = self.confDict['CATDB'] + '/tax'
        out_prefix = f'{wkdir}/{self.name}.CAT'
        orf_annot = f'{out_prefix}.ORF2LCA.txt'
        orf_rate = f'{out_prefix}.host_orf_rate.tsv'
        cmd.extend(
            ['CAT_pack contigs', '-c', in_fa, '-t', self.threads, '-d', CAT_DB, '-t', CAT_TAX, '-o', out_prefix, '\n',
            'host_orf_rate_from_CAT.py', orf_annot,orf_rate, '\n']
        )
        return cmd, wkdir
