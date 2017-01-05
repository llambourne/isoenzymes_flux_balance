import os
from collections import defaultdict
import pandas as pd


# TODO: process dnds data
# TODO: add FBA model downloads


def genetic_interaction_data(dataDir):
    download_and_filter_boone_data(
     'http://drygin.ccbr.utoronto.ca/~costanzo2009/sgadata_costanzo2009_intermediateCutoff_101120.txt.gz',
     '../data/processed/genetic_interactions.csv',
     dataDir)
    download_and_filter_boone_data(
     'http://drygin.ccbr.utoronto.ca/~costanzo2009/sgadata_costanzo2009_stringentCutoff_101120.txt.gz',
     '../data/processed/genetic_interactions_stringent.csv',
     dataDir)


def download_and_filter_boone_data(url, outpath, dataDir):
    zipFileName = url.split('/')[-1]
    fileName = zipFileName[:-3]
    if not os.path.exists(os.path.join(dataDir, fileName)):
        os.system('wget -P ' + dataDir + ' ' + url)
        os.system('gunzip ' + os.path.join(dataDir, zipFileName))
    # http://drygin.ccbr.utoronto.ca/~costanzo2009/
    columns = ['Query ORF',
               'Query gene name',
               'Array ORF',
               'Array gene name',
               'epsilon',
               'Standard deviation',
               'p-value']
    expData = pd.read_table('../data/external/'+fileName, names=columns)
    # drop rows with missing epsilon/p-values
    expData.dropna(axis=0, inplace=True)
    expData['interaction'] = (expData['epsilon'] > 0).map({True: 'positive',
                                                           False: 'negative'})
    expData = expData[['Query ORF', 'Array ORF', 'interaction']]

    if os.path.exists(outpath):
        print outpath, 'already exists, skipping'
    else:
        expData.to_csv(outpath, index=False)


def full_genetic_interaction_data(rawDir='../data/external'):
    url = 'http://drygin.ccbr.utoronto.ca/~costanzo2009/sgadata_costanzo2009_rawdata_101120.txt.gz'
    if not os.path.exists(os.path.join(rawDir, url.split('/')[-1][:-3])):
        os.system('wget -P ' + os.path.join(rawDir, url))
        os.system('gunzip ../data/external/sgadata_costanzo2009_rawdata_101120.txt.gz')
    # http://drygin.ccbr.utoronto.ca/~costanzo2009/
    columns = ['Query ORF',
               'Query gene name',
               'Array ORF',
               'Array gene name',
               'epsilon',
               'Standard deviation',
               'p-value',
               'Query SMF',
               'Query SMF std',
               'Array SMF',
               'Array SMF std',
               'Double mutant fitness',
               'Double mutant fitness std']
    expData = pd.read_table('../data/external/sgadata_costanzo2009_rawdata_101120.txt',
                            names=columns)
    # drop rows with missing epsilon/p-values
    expData.dropna(axis=0, inplace=True)

    def genetic_interaction(row):
        """Genetic interaction as defined in the paper.

        Args:
            row: row of a dataframe.

        Returns:
            str: positive/negative/none

        """
        if row['epsilon'] > 0.08 and row['p-value'] < 0.05:
            return 'positive'
        elif row['epsilon'] < -0.08 and row['p-value'] < 0.05:
            return 'negative'
        else:
            return 'none'

    expData['interaction'] = expData.apply(genetic_interaction, axis=1)

    print expData.shape, 'before'   # test
    # filter to remove gene knockouts that have better growth than wildtype
    expData = expData[((expData['Query SMF'] - 2 * expData['Query SMF std']) <= 1.0) &
                      ((expData['Array SMF'] - 2 * expData['Array SMF std']) <= 1.0)]
    print expData.shape, 'after'   # test


    expData = expData[['Query ORF', 'Array ORF', 'interaction']]
    outpath = '../data/processed/genetic_interactions_filtered.csv'
    expData.to_csv(outpath, index=False)


def essentiality_data(rawDir, outDir):
    """Download lists of essential yeast genes.

    Args:
        rawDir (str): directory to save external datasets.

    """
    url = 'http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt'
    fileName = url.split('/')[-1]
    if not os.path.exists(os.path.join(rawDir, fileName)):
        os.system('wget -P ' + rawDir + ' ' + url)
    with open(os.path.join(rawDir, fileName), 'r') as fIn:
        with open(os.path.join(outDir, fileName), 'w') as fOut:
            for l in fIn.readlines()[2:-4]:
                    fOut.write(l.split('\t')[1] + '\n')


def download_yeast_models(rawDir):
    # TODO: WRITE THIS FUNCTION
    return


def combine_dnds_data(inDir, outDir):
    """Average rank of evolutionary rate.

    Reads in dN/dS data from files in a directory called dnds
    calculated for s. cer. genes using different yeast species.
    """
    outPath = os.path.join(outDir, 'dnds_rank_scerevisiae.tsv')
    dnds = defaultdict(dict)
    dfs = {}
    for fileName in os.listdir(inDir):
        if (fileName.startswith('dnds_scerevisiae_') and
            fileName.endswith('.tsv')):
            ortholog = fileName[17:-4]
            dfs[ortholog] = pd.read_table(os.path.join(inDir, fileName), header=0)
            dfs[ortholog]['dnds_rank'] = dfs[ortholog]['dN/dS*'].rank()
    df = pd.concat(dfs.values())
    if sum(pd.isnull(df['Saccharomyces cerevisiae ORF'])) > 0 or sum(pd.isnull(df['dN/dS*'])) > 0:
        raise UserWarning('ERROR in agregating dnds data.')
    avrgRank = df.groupby('Saccharomyces cerevisiae ORF').mean()['dnds_rank'].rank().copy()
    avrgRank.sort_values(inplace=True)
    avrgRank.to_csv(outPath, sep='\t')


def main():
    rawDir = '../data/external'
    processedDir = '../data/processed'
    download_yeast_models(rawDir)
    essentiality_data(rawDir, processedDir)
    genetic_interaction_data(rawDir)
    full_genetic_interaction_data()
    combine_dnds_data(os.path.join(rawDir, 'dnds'), processedDir)


if __name__ == '__main__':
    main()
