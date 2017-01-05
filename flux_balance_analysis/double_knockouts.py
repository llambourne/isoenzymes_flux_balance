#!/share/pkg/python/2.7.7/install/bin/python

import sys
import os
import cPickle as pickle
import numpy as np
import cobra
from cobra.flux_analysis import double_gene_deletion
from fba_utils import *


def double_function_knockouts(model, gene, otherGenes):
    """

    Args:
        model (cobra.model): FBA model.
        gene (str): orf ID of gene to knockout.
        otherGenes list(str): list of other ORFs to knockout
                              together with gene.

    Returns:
        dict: growth values of double knockouts.

    """
    growthValues = {}
    for geneB in otherGenes:
        growthValues[(gene, geneB)] = double_function_deletion(model, gene, geneB)
    return growthValues


def double_function_deletion(model, geneA, geneB):
    """Knock out all reactions associated with two genes.

    Args:
        model (cobra.model): FBA model.
        geneA (str): ORF ID.
        geneB (str): ORF ID.

    Returns:
        float: growth value after double knockout.

    """
    if geneA == geneB:
        raise UserWarning('Genes must not be the same.')
    modelGenes = set([gene.id for gene in model.genes])
    if geneA not in modelGenes:
        raise UserWarning('ERROR: '+geneA+' not in model')
    if geneB not in modelGenes:
        raise UserWarning('ERROR: '+geneB+' not in model')
    toKO = [r for r in model.reactions if geneA in r.gene_reaction_rule
                                          or geneB in r.gene_reaction_rule]
    oldBounds = [(r.lower_bound, r.upper_bound) for r in toKO]
    for r in toKO:
        r.lower_bound = 0.
        r.upper_bound = 0.
    growth = model.optimize().f
    # restore model to original state
    for r, bounds in zip(toKO, oldBounds):
        r.lower_bound = bounds[0]
        r.upper_bound = bounds[1]
    return growth


def main():
    if len(sys.argv) != 3:
        raise UserWarning('Wrong number of arguments')
    index = int(sys.argv[1])
    modelPath = sys.argv[2]
    model = load_sd_minus_his(modelPath)
    modelName = modelPath.split('/')[-1][:-4]
    orderedGeneNames = sorted([g.id for g in model.genes])
    glcGrowthVals = double_gene_deletion(model, gene_list1=[orderedGeneNames[index]],
                                         gene_list2=orderedGeneNames[index+1:])
    flcGrowthVals = double_function_knockouts(model, orderedGeneNames[index],
                                              orderedGeneNames[index+1:])
    outDir = '../models/' + modelName + '/cluster_output'
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    glcOutPath = os.path.join(outDir, 'growth_double_gene_knockouts_' +
                                       str(index)+'.pkl')
    flcOutPath = os.path.join(outDir, 'growth_double_function_knockouts_' +
                                      str(index)+'.pkl')
    with open(glcOutPath, 'w') as f:
        pickle.dump(glcGrowthVals, f)
    with open(flcOutPath, 'w') as f:
        pickle.dump(flcGrowthVals, f)


if __name__ == '__main__':
    main()
