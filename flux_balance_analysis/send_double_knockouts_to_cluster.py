import os
from os import path
import cobra


def send_jobs(modelPath):
    modelName = modelPath.split('/')[-1][:-4]
    outDir = '../models/' + modelName + '/cluster_output'
    if not path.exists(outDir):
        os.makedirs(outDir)
    model = cobra.io.read_sbml_model(modelPath)
    orfs = [g.id for g in model.genes]
    cmd = 'qsub -V -b y python double_knockouts.py '
    for i in range(len(model.genes) - 1):
        if not path.exists(path.join(outDir, 'growth_double_gene_knockouts_' 
                                             +str(i) + '.pkl')):
            os.system(cmd + ' ' + str(i) + ' ' + modelPath)


def main():
    send_jobs('../data/external/yeast_7.6/yeast_7.6.xml')


if __name__ == '__main__':
    main()
