import sys
from ..general import utils
from ..data.bioseq import Seq
from .vCont import EnviComp

class VirTaxa(Seq):
    '''
    Taxonomic assignment
    '''
    def __init__(self,fasta='',outdir='',threads=8):
        super().__init__(fasta,outdir)
        self.threads=str(threads)
    def ncbiRefSeqTaxa(self,orf_f):
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
        sed_cmd="'1i\QueryID\\tNCBI_ID\\tIdentity\\tAlnLen\\tMismatches\\tGap\\tQstart\\tQend\\tSstart\\tSend\\tEValue\\tBitScore'"
        cmd=[utils.selectENV('VC-General')]
        cmd.extend(
            ['blastp','-query',orf_f,'-out',blast_resu,
            '-db',refdb,'-num_threads',self.threads,
            '-outfmt 6 -evalue 1e-3 -max_target_seqs 1\n',
            "awk '$12>=50'",blast_resu,'>',filt_blast,'\n',
            'sed',sed_cmd,filt_blast,'>',filt_h_blast,'\n',
            "csvtk join -t -f 'NCBI_ID,NCBI_ID'",
            filt_h_blast,taxadb,'>',gene_taxa_blast,'\n',
            'virus_tax.py',orf_f,gene_taxa_blast,'>',votu_taxa,'\n']
        )
        return cmd,votu_taxa
    def demovir(self):
        '''
        Classify the virus contig by Demovir software for a certain single group.
        '''
        wkdir=f'{self.outdir}/demovir'
        utils.mkdir(wkdir)
        cmd=[utils.selectENV('VC-General')]
        cmd.extend(['cd',wkdir,'\n'])
        demovir_db=self.confDict['DemovirDB']
        TrEMBL_viral_taxa=f'{demovir_db}/TrEMBL_viral_taxa.RDS'
        cmd.extend(['ln -s',TrEMBL_viral_taxa,'\n'])
        uniprot_trembl_viral=f'{demovir_db}/uniprot_trembl.viral.udb'
        cmd.extend(['ln -s',uniprot_trembl_viral,'\n'])
        demovir=f'{sys.path[0]}/bin/demovir.*'
        cmd.extend(
            ['cp',demovir,'.\n',
            './demovir.sh',self.fasta,self.threads,'\n']
        )
        votu_taxa=f'{wkdir}/DemoVir_assignments.txt'
        return cmd,votu_taxa
    def mergeTaxa(self,taxa1,taxa2):
        taxa=f'{self.outdir}/{self.name}.votu.taxa.txt'
        cmd=[utils.selectENV('VC-General')]
        cmd.extend(['merge_taxa.pl',taxa1,taxa2,'>',taxa,'\n'])
        return cmd
    def Classify(self,unrun=False):
        cmd=[self.envs]
        tmp_cmd,orf_faa=self.genePred()
        cmd.extend(tmp_cmd)
        tmp_cmd,ncbi_taxa=self.ncbiRefSeqTaxa(orf_faa)
        cmd.extend(tmp_cmd)
        tmp_cmd,demovir_taxa=self.demovir()
        cmd.extend(tmp_cmd)
        Comp=EnviComp(
            fasta=orf_faa,
            outdir=self.outdir,
            threads=self.threads
        )
        cmd.extend(Comp.vContact(''))
        cmd.extend(self.mergeTaxa(demovir_taxa,ncbi_taxa))
        shell=f'{self.outdir}/{self.name}_classify.sh'
        utils.printSH(shell,cmd)
        results=''
        if not unrun: results=utils.execute(shell)
        return results
