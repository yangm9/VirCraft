from os import path
from ..process import cmdExec, general
from ..identify.viridsop import VirSurvey

class VirRef(VirSurvey):
    '''

    '''
    def __init__(self, config, outdir):
        VirSurvey.__init__(self, config, outdir)
        self.datadir = self.wkdir #Save the work directory for last step
        self.wkdir = f'{self.outdir}/03.vOTUs' #Update the Work directory
        self.votus = f'{self.wkdir}/all_votus.fa'
        self.qual_summ = f'{self.wkdir}/2.checkv/quality_summary.tsv' #Used for abundance-vs-length scatter plot in VirTPM.
        general.mkdir(self.wkdir)
    def seqCluster(self, fasta: str, group: str):
        '''
        Cluster the sequence and remove redundancy for FastA file.
        '''
        cdhit_cmd = [self.envs]
        votus_fa = f'{self.wkdir}/{group}_votus.fa'
        cdhit_cmd.extend(
            ['cd-hit-est', '-i', fasta, '-o', votus_fa, 
            '-c 0.95 -aS 0.85 -n 10 -d 0 -M 160000 -T 28\n']
        )
        checkv_dir = f'{self.wkdir}/2.checkv'
        general.mkdir(checkv_dir)
        checkv_db = self.confDict["CheckVDB"]
        cdhit_cmd.extend(
            ['checkv end_to_end', votus_fa, checkv_dir, '-d', checkv_db, '-t 40\n']
        )
        cdhit_sh = f'{self.wkdir}/{group}_cdhit.sh'
        general.printSH(cdhit_sh, cdhit_cmd)
        results = cmdExec.execute(cdhit_cmd)
        return results
    def RmDup(self):
        results = ''
        merge_fa_cmd = ['cat']
        for grp in groups:
            fasta = f'{self.datadir}/{grp}/curation/virus_positive.fna'
            merge_fa_cmd.append(fasta)
            grpdir = f'{self.wkdir}/{grp}'
            general.mkdir(grpdir)
            results += f'{grp}: \n'
            results += seqCluster(fasta, grp, wkdir)
        all_fasta = f'{self.wkdir}/all_virus_positive.fa'
        merge_fa_cmd.append(f'>{all_fasta}\n')
        merge_fa_sh = f'{self.wkdir}/merge_fa.sh'
        general.printSH(merge_fa_sh, merge_fa_cmd)
        results += cmdExec.execute(merge_fa_cmd)
        results += self.seqCluster(merged_fasta, 'all', wkdir)
        return results
