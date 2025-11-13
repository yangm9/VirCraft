import os
from ..general import utils
from ..data.bioseq import VirSeq

class VirDetectTools(VirSeq):
    '''
    Generate command for VirSorter2, DeepVirFinder, VIBRANT, geNomad, Contig Annotation Tool (CAT).
    '''
    vs2_subcmds = ['--keep-original-seq', '--seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv']
    def __init__(self, fasta=None, outdir=None, threads=8):
        super().__init__(fasta, outdir)
        self.threads = str(threads)
    def virsorter(self, in_fa: str, n: int, min_len=1500, min_score=0.5): #n is the index corespoding 2 purpose: 0: for viral identification, and 1 for DRAMv input files
        idx = str(n+1)
        min_score = str(min_score)
        min_len = str(min_len)
        wkdir = f'{self.wkfile_dir}/vs2-pass{idx}'
        utils.mkdir(wkdir)
        cmd = [utils.selectENV('VC-VirSorter2')]
        cmd.extend(
            ['virsorter run', self.vs2_subcmds[n], '-i', in_fa, '-d', self.confDict['VirSorter2DB'], '-w', wkdir, '-j', self.threads, '--min-length', min_len, '--min-score', min_score, self.confDict['VirSorter2Opts'] + '\n']
        )
        return cmd, wkdir

    def deepvirfinder(self, min_len=1500):
        min_len = str(min_len)
        input_prefix = f'{self.wkfile_dir}/{self.name}'
        input_fasta = input_prefix + '.lt2100000.fa'
        cmd = self.lenCutoff(0, 2100000)
        cmd.extend([utils.selectENV('VC-DeepVirFinder')])
        wkdir = f'{self.wkfile_dir}/deepvirfinder'
        utils.mkdir(wkdir)
        cmd.extend(
            ['dvf.py', '-i', input_fasta, '-m', self.confDict['DeepVirFinderDB'], '-o', wkdir, '-c', self.threads, '-l', min_len, self.confDict['DeepVirFinderOpts'] + '\n']
        )
        return cmd, wkdir

    def vibrant(self, min_len=1500):
        min_len = str(min_len)
        cmd = [utils.selectENV('VC-VIBRANT')]
        wkdir = f'{self.wkfile_dir}/VIBRANT_{self.name}'
        VIBRANT_DB_databases = self.confDict['VIBRANTDB'] + '/databases'
        VIBRANT_DB_files = self.confDict['VIBRANTDB'] + '/files'
        cmd.extend(
            ['VIBRANT_run.py', '-i', self.fasta, '-l', min_len, '-t', self.threads, '-folder', self.wkfile_dir, '-d', VIBRANT_DB_databases, '-m', VIBRANT_DB_files, self.confDict['VIBRANTOpts'] + '\n']
        )
        return cmd, wkdir
    def genomad(self):
        cmd = [utils.selectENV('VC-geNomad')]
        wkdir = f'{self.wkfile_dir}/genomad'
        geNomad_DB = self.confDict['geNomadDB']
        cmd.extend(
            ['genomad end-to-end', '--threads', self.threads, self.fasta, wkdir, geNomad_DB, self.confDict['geNomadOpts'] + '\n']
        )
        return cmd, wkdir
    
    def contig_annotation_tool(self, in_fa=None):
        cmd = [utils.selectENV('VC-General')]
        wkdir = f'{self.wkfile_dir}/CAT_{self.name}'
        utils.mkdir(wkdir)
        CAT_DB = self.confDict['CATDB'] + '/db'
        CAT_TAX = self.confDict['CATDB'] + '/tax'
        out_prefix = f'{wkdir}/{self.name}.CAT'
        orf_annot = f'{out_prefix}.ORF2LCA.txt'
        orf_rate = f'{out_prefix}.host_orf_rate.tsv'
        cmd.extend(
            ['CAT_pack contigs', '-c', in_fa, '-t', self.threads, '-d', CAT_DB, '-t', CAT_TAX, '-o', out_prefix, self.confDict['CATOpts'] + '\n',
            'host_orf_rate_from_CAT.py', orf_annot,orf_rate, '\n']
        )
        return cmd, wkdir
