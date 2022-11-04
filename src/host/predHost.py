import sys
from ..process import cmdExec, general
from ..votus.deRep import VirRef

class VirHost():
    '''
    
    '''
    def __init__(self, config, outdir):
        VirRef.__init__(self, config, outdir)
        self.datadir = self.wkdir
        self.wkdir = f'{self.outdir}/08.hosts'
        general.mkdir(self.wkdir)
    def magsClassify(self):
        '''
        Classify the host MAGs by GTDBTK tools.
        '''
        host_mags = self.confDict['HostMAGs']
        cmd = [self.envs]
        cmd.extend(['ln -s', TrEMBL_viral_taxa, TrEMBL_viral_taxa_lnk, '\n'])
        uniprot_trembl_viral = f'{demovir_db}/uniprot_trembl.viral.udb'
        uniprot_trembl_viral_lnk = f'{self.wkdir}/uniprot_trembl.viral.udb'
        demovir_cmd.extend(['ln -s', uniprot_trembl_viral, uniprot_trembl_viral_lnk, '\n'])
        demovir = f'{sys.path[0]}/bin/demovir.*'
        demovir_cmd.extend(['cp', demovir, self.wkdir, '\n'])
        votus = f'{self.outdir}/03.vOTUs/merged_virus_positive_nodup.fa'
        demovir = f'{self.wkdir}/demovir.sh'
        demovir_cmd.extend([demovir, votus, '32'])
        demovir_sh = f'{self.wkdir}/classify_by_demovir.sh'
        general.printSH(demovir_sh, demovir_cmd)
        results = cmdExec.execute(demovir_cmd)
        return results
