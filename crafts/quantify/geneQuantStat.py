from ..general import utils
from .multiQuant import multiGeneCount

class GeneAbdStat(multiGeneCount):
    '''
    Gene abundance main class
    '''
    def __init__(self,samp_info='',fasta='',outdir='',threads=8):
        super().__init__(samp_info,fasta,outdir,threads)
    def mergeAbd(self):
        abd=f'{self.outdir}/all_merged_gene.sf'
        cmd=['merge_tpms.pl',self.samp_info,self.outdir,'tpm Gene\n']
        return cmd,abd
    def QuantStat(self):
        self.geneCountBySamp()
        cmd=[self.envs]
        multicmd = f'''\
cd {self.outdir}
sh {self.name}_salmonidx.sh
processes=()
ls *_gene_count.sh | while read i
do
     echo -e "nohup sh $i 1>$i.o 2>$i.e&" | bash
     processes+=($!)
done
for process in "${{processes[@]}}"
do
     wait $process
done
'''
        cmd.append(multicmd)
        tmp_cmd,abd=self.mergeAbd()
        cmd.extend(tmp_cmd)
        cmd.append('rm -rf *_gene_count.sh*\n')
        shell=f'{self.outdir}/{self.name}_gene_count.sh'
        utils.printSH(shell,cmd)
        #results=utils.execute(cmd)
        return 0#results
