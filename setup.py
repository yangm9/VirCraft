import setuptools

setuptools.setup(name='VirCraft',
                 version='0.0.6',
                 description='VirCraft',
                 url='https://github.com/yangm9/ncbiDigger',
                 author=['Ming Yang'],
                 author_email='yangm@idsse.ac.cn',
                 license='GPLv3',
                 packages=setuptools.find_packages(),
                 package_data={'data': ['tRNAs/00_Promiscuous_tRNAs.fasta',
                                        'tRNAs/00_Promiscuous_tRNAs.nhr',
                                        'tRNAs/00_Promiscuous_tRNAs.nin',
                                        'tRNAs/00_Promiscuous_tRNAs.nsq',
                                        'decoy_viruses/Decoy_viruses_host_tax.tsv',
                                        'decoy_viruses/Decoy_viral_genomes.tar.gz'
                                        ]},
                 scripts=['bin/abd_by_taxa.py',
                          'bin/alpha_diversity.R',
                          'bin/barplot_for_taxa_tpm.R',
                          'bin/blast_classify.py'],
                 install_requires=[
                     'biopython>=1.72',
                     'pandas>=1.1.3',
                     'trnascan-se>=2.0.7',
                     'diamond==2.0.14.152',
                     'minced>=0.4.2',
                     'r-here>=1.0.0',
                     'r-seqinr>=4.2_5',
                     'r-dplyr>1.0.1',
                     'r-data.table>=1.13',
                     'r-stringr>=1.4.0',
                     'psutil>=5.7.2',
                     'r-lifecycle>=1.0.0'
                 ]
                 )
