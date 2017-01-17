import sys
import os

"""
Sends jobs to the computer cluster that calculate the cost for every gene
knockout for a simulation of a particular minimal media made up of one
metabolite from the carbon source list and one from the nitrogen source list.
"""


def main():
    modelPath = '../data/external/yeast_7.6/yeast_7.6.xml'
    with open('../data/processed/carbon_sources.txt', 'r') as f:
        carbonSourceNames = [l.strip() for l in f.readlines()]
    with open('../data/processed/nitrogen_sources.txt', 'r') as f:
        nitrogenSourceNames = [l.strip() for l in f.readlines()]
    cmd = 'qsub -V -b y "python single_knockouts.py '
    outDir = '../models/' + modelPath.split('/')[-1][:-4] + '/cluster_output'
    for carbon in carbonSourceNames:
        for nitrogen in nitrogenSourceNames:
            outPath = os.path.join(outDir, ('gene_loss_cost_' +
                            carbon + '_AND_' + nitrogen + '.pkl').replace(' ', '_'))
            args = modelPath + ' \\"' + carbon + '\\" \\"' + nitrogen + '\\""'
            if not os.path.exists(outPath):
                os.system(cmd + args)


if __name__ == '__main__':
    main()
