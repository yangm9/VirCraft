import sys
from optparse import OptionParser #该模块已经不在开发维护，改为argparse
from . import readme

class Options:
    def __init__(self,name:str,version:str):
        self.usage=readme.description(name,version)
        self.parser=OptionParser(self.usage,version=f'%prog {version}')
        self.parser.add_option(
            '-t','--threads', action='store', type='str',
            dest='threads',metavar='INT', default=False,
            help='Threads number.'
        )
        self.parser.add_option(
            '-o','--outdir', action='store', type='str',
            dest='outdir',metavar='STR', default=False,
            help='Output direcortory.'
        )
        if name=='assembly':
            self.parser.add_option(
                '-1','--fastq1', action='store', type='str',
                dest='fq1',metavar='STR', default=False,
                help='FastQ file for read 1.'
            )
            self.parser.add_option(
                '-2','--fastq2', action='store', type='str',
                dest='fq2',metavar='STR', default=False,
                help='FastQ file for read 2.'
            )
        elif name=='':
            pass

        self.opts,self.args=self.parser.parse_args()
    def checkOpts(self):
        if self.opts.outdir:
            outdir=os.path.abspath(self.opts.outdir)
            general.mkdir(outdir)
        else:
            self.parser.print_help()
            exit(0)
        return 0
