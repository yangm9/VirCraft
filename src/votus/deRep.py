from os import path
from ..general import cmdExec,general
from ..identify.viridsop import VirSurvey

class VirRef(VirSurvey):
    '''

    '''
    def __init__(self,config,outdir):
        VirSurvey.__init__(self,config,outdir)
        self.datadir=self.wkdir #Save the work directory for last step
        self.wkdir=f'{self.outdir}/03.vOTUs' #Update the Work directory
        general.mkdir(self.wkdir)
        self.votus=f'{self.wkdir}/all_votus.fa'
        self.qual_summ=f'{self.wkdir}/2.checkv/quality_summary.tsv' #Used for abundance-vs-length scatter plot in VirTPM.
    def seqCluster(self,fasta: str,group: str):
        '''
        Cluster the sequence and remove redundancy for FastA file.
        '''
        cmd=[self.envs]
        votus_fa=f'{self.wkdir}/{group}_votus.fa'
        cmd.extend(
            ['cd-hit-est','-i',fasta,'-o',votus_fa,
            '-c 0.95 -aS 0.85 -n 10 -d 0 -M 160000 -T 28\n']
        )
        wkdir=f'{self.wkdir}/2.checkv'
        general.mkdir(wkdir)
        checkv_db=self.confDict['CheckVDB']
        cmd.extend(
            ['checkv end_to_end',votus_fa,wkdir,'-d',checkv_db,'-t 40\n']
        )
        shell=f'{self.wkdir}/{group}_cdhit.sh'
        general.printSH(shell,cmd)
        results=cmdExec.execute(cmd)
        return results
    def RmDup(self):
        results=''
        cmd=['cat']
        for grp in groups:
            fasta=f'{self.datadir}/{grp}/curation/virus_positive.fna'
            merge_fa_cmd.append(fasta)
            wkdir=f'{self.wkdir}/{grp}'
            general.mkdir(grpdir)
            results+=f'{grp}: \n'
            results+=seqCluster(fasta,grp,wkdir)
        all_fasta=f'{self.wkdir}/all_virus_positive.fa'
        cmd.append(f'>{all_fasta}\n')
        shell=f'{self.wkdir}/merge_fa.sh'
        general.printSH(shell,cmd)
        results+=cmdExec.execute(merge_fa_cmd)
        results+=self.seqCluster(merged_fasta,'all',wkdir)
        return results
