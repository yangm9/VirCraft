from ..general import utils
from ..data.bioseq import Seq

class GeneFunc(Seq):
    '''
    Gene function annotation.
    1) Predict genes from fasta file;
    2) Annotate the function of these genes.
    '''
    def __init__(self, fasta=None, outdir=None, threads=8):
        super().__init__(fasta, outdir)
        self.threads = str(threads)
    def eggnogAnno(self, orfs_faa):
        wkdir = f'{self.outdir}/eggnog'
        utils.mkdir(wkdir)
        anno_prefix = f'{wkdir}/{self.name}'
        seed_orth = f'{anno_prefix}.emapper.seed_orthologs'
        eggout = f'{wkdir}/{self.name}_eggout'
        eggout_tsv = f'{eggout}/.emapper.annotations'
        eggnog_db = self.confDict['EggNOGDB']
        cmd = [utils.selectENV('VC-General')]
        cmd.extend(
            ['emapper.py -m diamond --no_annot --no_file_comments', '--cpu', self.threads, '-i', orfs_faa, '-o', anno_prefix, '--data_dir', eggnog_db, '\n', 
             'emapper.py', '--annotate_hits_table', seed_orth, '--no_file_comments', '-o', eggout, '--cpu', self.threads, '--data_dir', eggnog_db, '--override\n', 
             'cd', wkdir, '&& GO_anno_from_tab.py', '-i', eggout_tsv, '&& cd -\n']
        )
        return cmd
    def keggAnno(self, orfs_faa):
        wkdir = f'{self.outdir}/kegg'
        utils.mkdir(wkdir)
        kegg_db = self.confDict['KofamscanDB']
        ko_prof = f'{kegg_db}/profiles'
        ko_list = f'{kegg_db}/ko_list'
        exec_anno = f'{wkdir}/{self.name}.exec_annotation.tsv'
        exec_anno_detail = f'{wkdir}/{self.name}.exec_annotation.detail.tsv'
        kegg_info = f'{wkdir}/{self.name}.exec_annotation.kegg_anno.tsv'
        kegg_lv2_stat = f'{wkdir}/kegg_lv2_stat.tsv'
        cmd = [utils.selectENV('VC-General')]
        cmd.extend(
            ['exec_annotation -f mapper', '--cpu', self.threads, '-p', ko_prof, '-k', ko_list, '-o', exec_anno, orfs_faa, '\n', 
             'exec_annotation -f detail', '--cpu', self.threads, '-p', ko_prof, '-k', ko_list, '-o', exec_anno_detail, orfs_faa, '\n', 
             'cd', wkdir, '&& kaas_kofam2pathwayAnalysis.py', exec_anno, '&& cd -\n', 
             'kegg_lv_stat.py', kegg_info, 'level_2 >', kegg_lv2_stat, '\n', 
             'kegg_lv2_stat.R', kegg_lv2_stat, wkdir, '\n']
        )
        return cmd
    def FuncAnnot(self):
        tmp_cmd, orfs_faa = self.genePred()
        cmd.extend(tmp_cmd)
        cmd.extend(self.eggnogAnno(orfs_faa))
        cmd.append(utils.selectENV('VC-General'))
        cmd.extend(self.keggAnno(orfs_faa))
        shell = f'{self.outdir}/{self.name}_gene_annot.sh'
        utils.printSH(shell, cmd)
        results = utils.execute(shell)
        return results
