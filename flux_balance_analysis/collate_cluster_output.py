"""Script to read in all the costs and write to csv."""

import os
import cPickle as pickle
from collections import defaultdict
import cobra
from fba_utils import *


def collate_cluster_output(modelPath):
    modelName = modelPath.split('/')[-1][:-4]
    modelDir = '../models/' + modelName
    clusterDir = os.path.join(modelDir, 'cluster_output')
    flcOutPath = os.path.join(modelDir, 'function_loss_costs.tsv')
    glcOutPath = os.path.join(modelDir, 'gene_loss_costs.tsv')
    noGrowthOutPath = os.path.join(modelDir, 'minimal_media_no_growth.tsv')
    dblGenePath = os.path.join(modelDir, 'double_gene_loss_growth.csv')
    dblFunctionPath = os.path.join(modelDir, 'double_function_loss_growth.csv')
    model = cobra.io.read_sbml_model(modelPath)
    modelGeneNames = set([g.id for g in model.genes])
    with open('../data/processed/carbon_sources.txt', 'r') as f:
        carbonSourceNames = [l.strip() for l in f.readlines()]
    with open('../data/processed/nitrogen_sources.txt', 'r') as f:
        nitrogenSourceNames = [l.strip() for l in f.readlines()]
    mediaWithGrowth = []
    mediaNoGrowth = []
    glc = defaultdict(list)
    flc = defaultdict(list)
    dblGLC = {}
    dblFLC = {}
    for c in carbonSourceNames:
        for n in nitrogenSourceNames:
            fnGLC = ('gene_loss_cost_' + c + '_AND_' + n + '.pkl').replace(' ', '_')
            glcInPath = os.path.join(clusterDir, fnGLC)
            if not os.path.exists(glcInPath):
                model = minimal_media(model, c, n)
                if model.optimize().f > 0.01:
                    raise UserWarning('growth on media but no file found\n' +
                                       c + ' AND ' + n)
                mediaNoGrowth.append((c, n))
                continue
            with open(glcInPath, 'r') as f:
                mediaWithGrowth.append((c, n))
                glcOneMedia = pickle.load(f)
                if set(glcOneMedia.keys()) != modelGeneNames:
                    raise UserWarning('Unknown or missing genes in:'+glcInPath)
                for geneName, cost in glcOneMedia.items():
                    glc[geneName].append(cost)
            fnFLC = ('function_loss_cost_' + c + '_AND_' + n + '.pkl').replace(' ', '_')
            flcInPath = os.path.join(clusterDir, fnFLC)
            with open(flcInPath, 'r') as f:
                flcOneMedia = pickle.load(f)
                if set(flcOneMedia.keys()) != modelGeneNames:
                    raise UserWarning('Unknown or missing genes in:'+flcInPath)
                for geneName, cost in flcOneMedia.items():
                    flc[geneName].append(cost)
    print len(mediaNoGrowth), '/', len(carbonSourceNames) * len(nitrogenSourceNames),
    print 'Environments with no wildtype growth'
    ########### Load double function knockouts #######################
    for i in range(len(model.genes) - 1):
        filePath = os.path.join(clusterDir,
                                'growth_double_function_knockouts_'+str(i)+'.pkl')
        if not os.path.exists(filePath):
            raise UserWarning(filePath + ' does not exist')
        with open(filePath, 'r') as f:
                    dblFunctionKO = pickle.load(f)
                    if len(dblFunctionKO) != (len(model.genes) - i) - 1:
                        raise UserWarning('Unexpected number of double knockouts'
                                          ' in file: ' + filePath + 'found ' +
                                           str(len(dblFunctionKO)) + ' expected: '
                                           + str((len(model.genes) - i) - 1))
                    if any([k in dblFLC for k in dblFunctionKO]):
                        raise UserWarning('Duplicate entries in results')
                    dblFLC.update(dblFunctionKO)
    ########### Load double gene knockouts #######################
    for i in range(len(model.genes) - 1):
        filePath = os.path.join(clusterDir,
                                'growth_double_gene_knockouts_'+str(i)+'.pkl')
        if not os.path.exists(filePath):
            raise UserWarning(filePath + ' does not exist')
        with open(filePath, 'r') as f:
            result = pickle.load(f)
            for i in range(result['data'].shape[0]):
                geneA = result['x'][i]
                for j in range(result['data'].shape[1]):
                    geneB = result['y'][j]
                    dblGLC[(geneA, geneB)] = result['data'][i, j]
    ################### write out data ###############################
    if not os.path.exists(noGrowthOutPath):
        with open(noGrowthOutPath, 'w') as f:
            for c, n in mediaNoGrowth:
                f.write(c + '\t' + n + '\n')
    else:
        print noGrowthOutPath, 'already exists'
    if not os.path.exists(flcOutPath):
        with open(flcOutPath, 'w') as f:
            f.write('ORF ID\t' + '\t'.join([' AND '.join(cn)
                                          for cn in mediaWithGrowth]) + '\n')
            for geneName, costs in flc.items():
                f.write(geneName + '\t' + '\t'.join([str(c) for c in costs])+'\n')
    else:
        print flcOutPath, 'already exists'
    if not os.path.exists(glcOutPath):
        with open(glcOutPath, 'w') as f:
            f.write('ORF ID\t' + '\t'.join([' AND '.join(cn)
                                          for cn in mediaWithGrowth]) + '\n')
            for geneName, costs in glc.items():
                f.write(geneName + '\t' + '\t'.join([str(c) for c in costs])+'\n')
    else:
        print glcOutPath, 'already exists'
    if not os.path.exists(dblGenePath):
        with open(dblGenePath, 'w') as f:
            for (geneA, geneB), growth in dblGLC.items():
                f.write(geneA + ',' + geneB + ',' + str(growth) + '\n')
    else:
        print dblGenePath, 'already exists'
    if not os.path.exists(dblFunctionPath):
        with open(dblFunctionPath, 'w') as f:
            for (geneA, geneB), growth in dblFLC.items():
                f.write(geneA + ',' + geneB + ',' + str(growth) + '\n')
    else:
        print dblFunctionPath, 'already exists'


def main():
    collate_cluster_output('../data/external/yeast_7.6/yeast_7.6.xml')


if __name__ == '__main__':
    main()
