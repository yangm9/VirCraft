from ..general import utils
from ..data.bioseq import Reads
from ..data.bioseq import Seq

class Assembly(Reads):
    '''
    This class assembles the clean reads using SPAdes and/or MEGAHIT assembles.
    '''
    # Set the environment variables for the assembly process using a utility function
    envs=utils.selectENV('VC-Assembly')
    def __init__(self, fq1=None, fq2=None, outdir=None, threads=8):
        '''
        Initialize the Assembly class.
        :param fq1: Path to the FASTQ1 file.
        :param fq2: Path to the FASTQ2 file.
        :param outdir: output directory. .
        :param threads: Number of threads.
        '''
        # Initialize from the parent class Reads
        super().__init__(fq1, fq2, outdir)
        self.threads = str(threads)
        self.methDict = {
            'm': self.megahit, # 'm' key refers to the MEGAHIT assembly method
            's': self.spades   # 's' key refers to the SPAdes assembly method
        }

    def spades(self, fastqs: list):
        '''
        Assemble the metagenome by SPAdes.
        :param fastqs: List of FASTQ files.
        :return: The command to run SPAdes and the resulting scaffold file.
        '''
        wkdir = f'{self.wkfile_dir}/spades'
        utils.mkdir(wkdir)
        # Prepare input parameters based on the number of FASTQ files
        if len(fastqs) == 1:
            input_para = f'-s {fastqs[0]}'
        elif len(fastqs) == 2:
            input_para = f'--pe1-1 {fastqs[0]} --pe1-2 {fastqs[1]}'
        elif len(fastqs)==3:
            input_para = f'--pe1-1 {fastqs[0]} --pe1-2 {fastqs[1]} -s {fastqs[2]}'
        else:
            pass
        cmd = ['spades.py', input_para, '-t', self.threads, '-o', wkdir, self.confDict['SPAdesOpts'], '\n']
        scaf=f'{wkdir}/scaffolds.fasta'
        return cmd, scaf

    def megahit(self, fastqs: list):
        '''
        Assemble metagenome by megahit.
        :param fastqs: List of FASTQ files.
        :return: The command to run MEGAHIT and the resulting scaffold file.
        '''
        input_para = ''
        wkdir = f'{self.wkfile_dir}/megahit'
        tmpdir = f'{self.wkfile_dir}/megahit.tmp'
        utils.mkdir(tmpdir)

        # Prepare input parameters based on the number of FASTQ files
        if len(fastqs) == 1:
            input_para = f'-r {fastqs[0]}'
        elif len(fastqs) == 2:
            input_para = f'-1 {fastqs[0]} -2 {fastqs[1]}'
        elif len(fastqs) == 3:
            input_para = f'-1 {fastqs[0]} -2 {fastqs[1]} -r {fastqs[2]}'
        else:
            pass

        # Create the command for running MEGAHIT
        cmd = ['megahit', input_para, '-o', wkdir, '-t', self.threads, self.confDict['MegahitOpts'], '--tmp-dir', tmpdir, '\n']
        scaf=f'{wkdir}/final.contigs.fa'
        return cmd, scaf

    def unmapReads(self, scaf: str):
        '''
        Align the FASTQ reads back to the assembled contigs.
        :param scaf: Path to the assembled scaffold file.
        :return: The command to map reads back to the scaffold and the paths to unmapped reads.
        '''
        wkdir = f'{self.wkfile_dir}/alignment'
        utils.mkdir(wkdir)

        # Define paths for BWA index and output files
        bwa_idx = f'{wkdir}/scaffoldsIDX'
        unused_bam = f'{wkdir}/unused_reads.bam'
        unused_fq_1 = f'{wkdir}/unused_reads_1.fq'
        unused_fq_2 = f'{wkdir}/unused_reads_2.fq'
        unused_fq_s = f'{wkdir}/unused_reads_s.fq'
        # Create the command to index the scaffolds, map the reads, and extract unmapped reads
        cmd = ['bwa index -a bwtsw', scaf, '-p', bwa_idx, '\n',
            'bwa mem', '-t', self.threads, bwa_idx, self.fastqs[0], self.fastqs[1], '|samtools view -bf 4 >', unused_bam, '\n',
            'samtools fastq -N', unused_bam, '-1', unused_fq_1, '-2', unused_fq_2, '-s', unused_fq_s, '\n']
        # List of unused FASTQ files (unmapped reads)
        unused_fqs = [unused_fq_1, unused_fq_2, unused_fq_s]
        return cmd, unused_fqs

    def mixAsse(self, fastqs, process='m'):
        '''
        Perform the assembly process based on the method defined in 'process'.
        :param fastqs: List of FASTQ files.
        :param process: Process to follow, default is MEGAHIT ('m').
        :return: Assembly command and final scaffold file path.
        '''
        cmd = []
        scafs = []
        tmp_cmd, tmp_scaf = self.methDict[process[0]](fastqs)
        cmd.extend(tmp_cmd)
        scafs.append(tmp_scaf)
        final_scaf = f'{self.outdir}/final_assembly.fasta'
        steps = len(process)

        # If two methods are specified, process the assembly in two steps
        if steps == 2:
            tmp_cmd, unused_fqs = self.unmapReads(tmp_scaf)
            cmd.extend(tmp_cmd)
            tmp_cmd, tmp_scaf = self.methDict[process[1]](unused_fqs) # Execute the second assembly method
            cmd.extend(tmp_cmd)
            scafs.append(tmp_scaf)
            cmd.extend(['cat', scafs[0], scafs[1], '>', final_scaf, '\n'])
        elif steps == 1:
            cmd.extend(['ln -s', tmp_scaf, final_scaf, '\n'])
        else: 
            pass
        return cmd, final_scaf
    
    def Assemble(self, process='ms', cutoff=1500, unrun=False, clear=False):
        '''
        Perform the entire assembly process and return the result.
        :param process: Assembly process (default is MEGAHIT).
        :param cutoff: Minimum contig length to retain.
        :param unrun: If True, the assembly will not run but the commands will be generated.
        :param clear: If True, intermediate files will be deleted after assembly.
        :return: Assembly results.
        '''
        
        cmd = [self.envs] # Store the environment variables command
        
        # Perform the assembly process and get commands and scaffold
        tmp_cmd, scaf = self.mixAsse(self.fastqs, process)
        cmd.extend(tmp_cmd)
        
        # Create a FASTA sequence object and append statistics commands
        FastA = Seq(scaf, self.outdir)
        cmd.extend(FastA.statFA())
        cmd.extend(FastA.lenCutoff())
        
        # Optionally clear intermediate alignment files
        if clear and len(process) == 2:
            alndir = f'{self.wkfile_dir}/alignment'
            scafIdx = f'{alndir}/scaffoldsIDX*'
            cmd.extend(['rm -f', scafIdx, '\n'])
        shell = f'{self.shell_dir}/{self.samp}_assembly.sh'
        utils.printSH(shell, cmd)
        
        # Execute the shell script if not in unrun mode
        results=''
        if not unrun: results = utils.execute(shell)
        return results
