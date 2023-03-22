import sys
from ..general import utils
from ..data.bioseq import Seq

class VirTaxa(Seq):
    '''
    Taxonomic assignment
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.threads=str(threads)
    def ncbiRefSeqTaxa(self,orfs_f):
        '''
        ORFs predicated from Prodigal (v2.6.3) were subjected to BLASTp (E-value of < 0.001, bitscore ≥ 50) against the NCBI viral RefSeq database (https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/).
        '''
        wkdir=f'{self.outdir}/blast'
        utils.mkdir(wkdir)
        dbdir=self.confDict['NCBIvRefProtDB']
        refdb=f'{dbdir}/viral.1.protein'
        taxadb=f'{dbdir}/NCBI_viral_full_taxnomomy.txt'
        blast_resu=f'{wkdir}/{self.name}.votu.blast'
        filt_blast=utils.insLable(blast_resu,'filt')
        filt_h_blast=utils.insLable(filt_blast,'h')
        gene_taxa_blast=utils.insLable(filt_h_blast,'sp')
        votu_taxa=f'{wkdir}/{self.name}.votu.taxa.txt'
        sed_cmd="'1i\QueryID\tNCBI_ID\tIdentity\tAlnLen\tMismatches\tGap\tQstart\tQend\tSstart\tSend\tEValue\tBitScore'"
        cmd=['blastp','-query',orfs_f,'-out',blast_resu,
            '-db',refdb,'-num_threads',self.threads,
            '-outfmt 6 -evalue 1e-3 -max_target_seqs 1\n',
            "awk '$12>=50'",blast_resu,'>',filt_blast,'\n',
            'sed',sed_cmd,filt_blast,'>',filt_h_blast,'\n',
            "csvtk join -t -f 'NCBI_ID,NCBI_ID'",
            filt_h_blast,taxadb,'>',gene_taxa_blast,'\n',
            'viruse_tax.py',gene_taxa_blast,'>',votu_taxa,'\n']
        return cmd,votu_taxa
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
        wkdir=f'{self.outdir}/demovir'
        utils.mkdir(wkdir)
        demovir=f'{wkdir}/demovir.sh'
        cmd.extend([demovir,self.fasta,self.threads,'\n'])
        votu_taxa=f'{wkdir}/DemoVir_assignments.txt'
        return cmd,votu_taxa
    def mergeTaxa(self,taxa1,taxa2):
        taxa=f'{self.outdir}/{self.name}.votu.taxa.txt'
        cmd=['merge_taxa.pl',taxa1,taxa2,'>',taxa,'\n']
        return cmd
    def Classify(self):
        cmd=[self.envs]
        tmp_cmd,orf_faa=self.genePred()
        cmd.extend(tmp_cmd)
        tmp_cmd,ncbi_taxa=self.ncbiRefSeqTaxa(orf_faa)
        cmd.extend(tmp_cmd)
        tmp_cmd,demovir_taxa=self.demovir()
        cmd.extend(tmp_cmd)
        cmd.extend(self.mergeTaxa(demovir_taxa,ncbi_taxa))
        shell=f'{self.outdir}/{self.name}_classify.sh'
        utils.printSH(shell,cmd)
        results=utils.execute(cmd)
        return results
