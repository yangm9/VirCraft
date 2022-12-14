import setuptools

setuptools.setup(name='VirMatcher',
                 version='0.3.3',
                 description='VirMatcher',
                 url='https://bitbucket.org/MAVERIClab/VirMatcher',
                 author=['Benjamin Bolduc', 'Ahmed Zayed'],
                 author_email='bolduc.10@osu.edu',
                 license='GPLv3',
                 packages=setuptools.find_packages(),
                 package_data={'data': ['tRNAs/00_Promiscuous_tRNAs.fasta',
                                        'tRNAs/00_Promiscuous_tRNAs.nhr',
                                        'tRNAs/00_Promiscuous_tRNAs.nin',
                                        'tRNAs/00_Promiscuous_tRNAs.nsq',
                                        'decoy_viruses/Decoy_viruses_host_tax.tsv',
                                        'decoy_viruses/Decoy_viral_genomes.tar.gz'
                                        ]},
                 scripts=['bin/VirMatcher',
                          'bin/Aggregator_combine_results.R',
                          'bin/ComputeNullParameters_R2.py',
                          'bin/ResultsAggregator.py'],
                 install_requires=[
                     'biopython>=1.72',
                     'pandas>=1.1.3',
                     'trnascan-se>=2.0.7',
                     'blast>=2.2.30',
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
