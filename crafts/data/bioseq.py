import os
from ..general import utils
from ..config.config import VirCfg

class Reads(VirCfg):
    '''
    FastQ processing class.
    '''
    envs = utils.selectENV('VC-General')
    postfixes = [
        '_1.fastq', '_1.fastq.gz', '_1.fq', '_1.fq.gz',
        '_R1.fastq', '_R1.fastq.gz', '_R1.fq', '_R1.fq.gz',
        '.R1.fastq', '.R1.fastq.gz', '.R1.fq', '.R1.fq.gz',
        '_1.clean.fq.gz'
    ]
    def __init__(self, fq1=None, fq2=None, outdir=None, *args, **kwargs):
        super().__init__()
        fq1 = os.path.abspath(fq1)
        fq2 = os.path.abspath(fq2)
        self.fastqs = [fq1, fq2]
        self.basename_fq1 = os.path.basename(self.fastqs[0])
        self.basename_fq2 = os.path.basename(self.fastqs[1])
        self.outdir = os.path.abspath(outdir)
        self.samp = self.getSampName
        utils.mkdir(self.outdir)
        self.shell_dir, self.wkfile_dir, self.stat_dir = utils.make_general_dir(self.outdir)
    @property
    def getSampName(self):
        samp = ''
        for post in self.postfixes:
            if self.basename_fq1.endswith(post):
                samp = self.basename_fq1.replace(post, '')
        return samp

class Seq(VirCfg):
    '''
    Fasta processing class.
    '''
    envs = utils.selectENV('VC-General')
    def __init__(self, fasta=None, outdir=None, *args, **kwargs):
        super().__init__()
        basename_fa = os.path.basename(fasta)
        self.name = os.path.splitext(basename_fa)[0]
        self.fasta = os.path.abspath(fasta)
        self.outdir = os.path.abspath(outdir)
        utils.mkdir(self.outdir)
        self.shell_dir, self.wkfile_dir, self.stat_dir = utils.make_general_dir(self.outdir)
    
    #Make bwa index for votus.
    def mkBwaIdx(self):
        cmd = [utils.selectENV('VC-Quantify')]
        bwa_idx = f'{self.wkfile_dir}/{self.name}.BWAIDX'
        cmd.extend(['bwa index -a bwtsw', self.fasta, '-p', bwa_idx, '\n'])
        return cmd, bwa_idx

    def lenCutoff(self, min_len=2000, max_len=None):
        cmd = [self.envs]
        filt_prefix = f'{self.outdir}/{self.name}'
        tmp_cmd = ['SeqLenCutoff.pl', self.fasta, filt_prefix, str(min_len)]
        if max_len is not None:
            tmp_cmd.append(str(max_len))

        cmd.extend(tmp_cmd)
        cmd.append('\n')
        return cmd

    def statFA(self):
        cmd = ['#Statistics and plotting of FASTA files\n', self.envs]
        size_dist = f'{self.stat_dir}/fasta_size_distribution.pdf'
        len_gc_stat = f'{self.stat_dir}/fasta_size_gc_stat.tsv'
        nl50_stat = f'{self.stat_dir}/fasta_nl50_stat.tsv'
        cmd.extend(
            ['fasta_size_distribution_plot.py', self.fasta, '-o', size_dist, '-s 2000 -g 10 -t "Sequence Size Distribution"\n',
             'fasta_size_gc.py', self.fasta, '>', len_gc_stat, '&& variables_scatter.R', len_gc_stat, 'Length~GC', self.stat_dir, '\n',
             'stat_NL50.pl', self.fasta, nl50_stat, '\n\n']
        )
        return cmd
    
    def genePred(self):
        cmd = [self.envs]
        wkdir = f'{self.wkfile_dir}/prodigal'
        utils.mkdir(wkdir)
        orf_ffn = f'{wkdir}/{self.name}.ffn'
        orf_faa = f'{wkdir}/{self.name}.faa'
        orf_gff = f'{wkdir}/{self.name}.gff'
        temp_faa = f'{wkdir}/temp.orf.faa'
        temp_ffn = f'{wkdir}/temp.orf.ffn'
        cmd = ['prodigal', '-i', self.fasta, '-d', temp_ffn, '-a', temp_faa, '-o', orf_gff, '-f gff -p meta -m -q\n',
               'cut -f 1 -d \" \"', temp_faa, '>', orf_faa, '&& cut -f 1 -d \" \"', temp_ffn, '>', orf_ffn, '\n',
               f'rm -f {wkdir}/temp.*\n']
        return cmd, orf_ffn, orf_faa

class VirSeq(Seq):
    def __init__(self, fasta=None, outdir=None, *args, **kwargs):
        super().__init__(fasta, outdir, *args, **kwargs)
    def checkv(self, in_fa: str):
        wkdir = f'{self.wkfile_dir}/checkv'
        utils.mkdir(wkdir)
        cmd = [utils.selectENV('VC-CheckV')]
        cmd.extend(
            ['checkv', 'end_to_end', in_fa, wkdir, '-d', self.confDict['CheckvDB'], '-t', self.threads, '\n']
        )
        provir_fna = f'{wkdir}/proviruses.fna'
        vir_fna = f'{wkdir}/viruses.fna'
        merged_fna = f'{wkdir}/combined.fna'
        cmd.extend(
            ['cat', provir_fna, vir_fna, '>', merged_fna, '\n']
        )
        return cmd, wkdir

class CDS(Seq):
    envs = utils.selectENV('VC-Quantify')
    def __init__(self, fasta=None, outdir=None, *args, **kwargs):
        super().__init__(fasta, outdir, *args, **kwargs)
    @property
    def mkSalmonIdx(self):
        cmd = [self.envs]
        wkdir = f'{self.wkfile_dir}/salmonidx'
        utils.mkdir(wkdir)
        idx = f'{wkdir}' # A directory
        cmd.extend(
            ['salmon index', '-p 8 -k 31', '-t', self.fasta, '-i', wkdir, '\n']
        )
        return cmd, idx

class ORF(Seq):
    envs = utils.selectENV('VC-General')
    def __init__(self, fasta=None, outdir=None, *args, **kwargs):
        super().__init__(fasta, outdir, *args, **kwargs)
    def eggnogAnno(self):
        wkdir = f'{self.outdir}/eggnog'
        utils.mkdir(wkdir)
        anno_prefix = f'{wkdir}/{self.name}'
        seed_orth = f'{anno_prefix}.emapper.seed_orthologs'
        eggout = f'{wkdir}/{self.name}_eggout'
        eggnog_db = self.confDict['EggNOGDB']
        cmd = ['emapper.py -m diamond --no_annot --no_file_comments', '--cpu', self.threads, '-i', self.fasta, '-o', anno_prefix, '--data_dir', eggnog_db, '\n',
               'emapper.py', '--annotate_hits_table', seed_orth, '--no_file_comments', '-o', eggout, '--cpu', self.threads, '--data_dir', eggnog_db, '--override\n']
        return cmd
