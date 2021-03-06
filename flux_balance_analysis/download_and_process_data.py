import os
from collections import defaultdict
from warnings import filterwarnings
import pandas as pd
import cobra
from single_knockouts import get_exchange_reactions

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
    outpath = '../data/processed/genetic_interactions_filtered.csv'
    if os.path.exists(outpath):
        print outpath, 'already exists. Skipping.'
        return
    url = 'http://drygin.ccbr.utoronto.ca/~costanzo2009/sgadata_costanzo2009_rawdata_101120.txt.gz'
    if not os.path.exists(os.path.join(rawDir, url.split('/')[-1][:-3])):
        os.system('wget -P ' + rawDir + ' ' + url)
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
    nBefore = expData.shape[0]
    # filter to remove gene knockouts that have better growth than wildtype
    expData = expData[((expData['Query SMF'] - 2 * expData['Query SMF std']) <= 1.0) &
                      ((expData['Array SMF'] - 2 * expData['Array SMF std']) <= 1.0)]
    print 'removed', nBefore - expData.shape[0], 'out of', nBefore
    print 'knockouts with higher growth than wildtype'
    expData = expData[['Query ORF', 'Array ORF', 'interaction']]
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


def download_yeast_model(rawDir):
    """Download FBA model."""
    url = 'https://pilotfiber.dl.sourceforge.net/project/yeast/yeast_7.6.zip'
    outPath = os.path.join(rawDir, 'yeast_7.6')
    if os.path.exists(outPath):
        print outPath, 'already exists. Skipping.'
        return
    os.system('wget -P ' + rawDir + ' ' + url)
    zipPath = os.path.join(rawDir, url.split('/')[-1])
    os.system('unzip ' + zipPath + ' -d ' + rawDir)
    os.remove(zipPath)


def get_blocked_reactions(modelPath):
    """Write out list of blocked reactions.

    A blocked reaction is defined as one for with no flux
    can go through in the case when all exchange reactions
    (i.e. all intake of nutrients) are open.

    Args:
        modelPath (str): path to FBA model xml file.

    """
    outDir = '../models/yeast_7.6'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    outPath = os.path.join(outDir, 'blocked_genes.txt')
    if os.path.exists(outPath):
        print outPath, 'already exists. Skipping.'
        return
    filterwarnings('ignore', 'charge of s_[0-9][0-9][0-9][0-9] is not a number ()')
    filterwarnings('ignore', 'uppercase AND/OR found in rule ')
    model = cobra.io.read_sbml_model(modelPath)
    exchangeReactions = get_exchange_reactions(model)
    for r in exchangeReactions:
        r.lower_bound = -10.
    blockedReactions = []
    for reaction in model.reactions:
        model.objective = reaction.id
        if model.optimize(objective_sense='maximize').f == 0. and \
           model.optimize(objective_sense='minimize').f == 0.:
            blockedReactions.append(reaction.id)
    print len(blockedReactions), 'blocked reactions'
    blockedGenes = []
    for gene in model.genes:
        if all([r.id in blockedReactions for r in gene.reactions]):
            blockedGenes.append(gene.id)
    print len(blockedGenes), 'blocked genes'
    with open(outPath, 'w') as f:
        f.write('\n'.join(blockedGenes))


def combine_dnds_data(inDir, outDir):
    """Average rank of evolutionary rate.

    Reads in dN/dS data from files in a directory called dnds
    calculated for s. cer. genes using different yeast species.
    """
    outPath = os.path.join(outDir, 'dnds_rank_scerevisiae.tsv')
    if os.path.exists(outPath):
        print outPath, 'already exists. Skipping.'
        return
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
    download_yeast_model(rawDir)
    essentiality_data(rawDir, processedDir)
    genetic_interaction_data(rawDir)
    full_genetic_interaction_data()
    combine_dnds_data(os.path.join(rawDir, 'dnds'), processedDir)
    get_blocked_reactions(os.path.join(rawDir, 'yeast_7.6/yeast_7.6.xml'))


if __name__ == '__main__':
    main()
