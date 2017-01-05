#!/share/pkg/python/2.7.7/install/bin/python

import sys
import os
import cPickle as pickle
import numpy as np
import cobra
from cobra.flux_analysis import single_gene_deletion, single_reaction_deletion
from fba_utils import *


def single_knockout_modified_loss_cost(model):
    """A modified gene-loss cost that assumes isoenzymes are non-redundant.

    Args:
        model (cobra.model): FBA model.

    Returns:
        float: growth value after double knockout.

    """
    wtGrowth = model.optimize(solver='gurobi').f
    lossCosts = {}
    for orfID in [g.id for g in model.genes]:
        toKO = [r for r in model.reactions if orfID in r.gene_reaction_rule]
        oldBounds = [(r.lower_bound, r.upper_bound) for r in toKO]
        for r in toKO:
            r.lower_bound = 0.
            r.upper_bound = 0.
        growth = model.optimize(solver='gurobi').f
        # restore model to original state
        for r, bounds in zip(toKO, oldBounds):
            r.lower_bound = bounds[0]
            r.upper_bound = bounds[1]
        lossCosts[orfID] = (wtGrowth - growth) / wtGrowth
    return lossCosts


def single_knockout_loss_costs(model, verbose=False):
    """

    Args:
        model (cobra.model): FBA model.
        verbose (bool): flag to print out messages.

    Returns:
        bool: did wildtype grow in these conditions.
        dict(str: list(float)): gene-loss cost per gene.
        dict(str: list(float)): function-loss cost per gene.

    """
    reactionsWithGenesMapped = [r for r in model.reactions if r.gene_reaction_rule != '']
    wtGrowth = model.optimize(solver='gurobi').f
    if wtGrowth < 0.01:
        if verbose:
            print 'wildtype failed to grow'
        return False, {}, {}
    if verbose:
        print 'running gene deletions'
    geneDelGrowth, geneDelStatus = single_gene_deletion(model, solver='gurobi')
    if verbose:
        print 'finished gene deletions, running reaction deletions'
    reactionDelGrowth, reactionDelStatus = single_reaction_deletion(model,
                                        reaction_list=reactionsWithGenesMapped,
                                        solver='gurobi')
    genes = geneDelGrowth.keys()
    geneLossCost = {gene: (wtGrowth - growth) / wtGrowth for gene, growth in geneDelGrowth.items()}
    ruleLossCost = {model.reactions.get_by_id(reaction).gene_reaction_rule: (wtGrowth - growth) / wtGrowth
                    for reaction, growth in reactionDelGrowth.items()}
    funcLossCost = {g: 0. for g in genes}
    for rule, cost in ruleLossCost.items():
        for gene in genes_in_rule(rule):
            funcLossCost[gene] += cost
    return True, geneLossCost, funcLossCost


def main():
    if len(sys.argv) != 4:
        raise UserWarning('ERROR: wrong number of arguments')
    modelPath = sys.argv[1]
    carbonSource = sys.argv[2]
    nitrogenSource = sys.argv[3]
    modelName = modelPath.split('/')[-1][:-4]
    model = minimal_media_model(modelPath, carbonSource, nitrogenSource)
    grew, glc, flc = single_knockout_loss_costs(model, verbose=False)
    if grew:
        outDir = '../models/' + modelName + '/cluster_output'
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        suffix = (carbonSource + '_AND_'
                  + nitrogenSource + '.pkl').replace(' ', '_')
        with open(os.path.join(outDir, 'gene_loss_cost_' + suffix), 'w') as f:
            pickle.dump(glc, f)
        with open(os.path.join(outDir, 'function_loss_cost_' + suffix), 'w') as f:
            pickle.dump(flc, f)


if __name__ == '__main__':
    main()
