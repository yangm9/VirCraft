#!/usr/bin/env python3
import sys
import pandas as pd

def sumAbundByTaxa(TaxaAbund):
    '''
    Sum the abundance by Taxa for each sample.
    '''
    df = pd.read_csv(TaxaAbund_f, sep='\t', header=0)
    df = df.iloc[:, 0:-2]
    TaxaList = df['Contig'].unique().tolist()
    df = df.set_index('Contig')
    sumDF = pd.DataFrame()
    for i in range(len(TaxaList)):
        series = df.loc[TaxaList[i]]
        if isinstance(series,pd.DataFrame):
            series = pd.Series(series.sum(),name=TaxaList[i])
        sumDF = sumDF.append(series)

    sumDF.reset_index(inplace=True)
    sumDF.rename(columns={'index':'Contig'}, inplace=True)
    sumDF.fillna(0, inplace=True)
    sumDF.to_csv('tax_tpm.xls', sep='\t',index=False)
    return 0

if __name__=='__main__':
    sumAbundByTaxa(sys.argv[1],sys.argv[2])
