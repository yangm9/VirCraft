import setuptools

setuptools.setup(name='VirCraft',
                 version='0.0.1',
                 description='VirCraft',
                 url='https://github.com/yangm9/VirCraft',
                 author=['Ming Yang'],
                 author_email='yangm@idsse.ac.cn',
                 license='GPLv3',
                 packages=setuptools.find_packages(),
                 package_data={'data': ['tRNAs/00_Promiscuous_tRNAs.fasta',
                                        ]},
                 scripts=['bin/VirMatcher',
                          'bin/Aggregator_combine_results.R',
                          ],
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
