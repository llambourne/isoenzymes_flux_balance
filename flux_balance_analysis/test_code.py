"""
Sanity checks on data and code.
"""

import os
import numpy as np
import cobra
from enzyme import enzyme
from warnings import filterwarnings

class TestFBAModel:

    def setup_class(self):
        modelPath='../data/external/yeast_7.6/yeast_7.6.xml'
        filterwarnings('ignore', 'charge of s_[0-9][0-9][0-9][0-9] is not a number ()')
        filterwarnings('ignore', 'uppercase AND/OR found in rule ')
        self.model = cobra.io.read_sbml_model(modelPath)
        modelDir = '../models/yeast_7.6'
        self.genes = {g.id: enzyme(g.id) for g in self.model.genes}
        for g in self.model.genes:
            self.genes[g.id].reactionRules = [r.gene_reaction_rule for r in g.reactions]
        with open(os.path.join(modelDir, 'gene_loss_costs.tsv'), 'r') as f:
            lines = f.readlines()
            minimalMedia = [tuple(m.split(' AND ')) for m in lines[0].strip().split('\t')[1:]]
            for line in lines[1:]:
                self.genes[line.split('\t')[0]].geneLossCosts = np.array([float(i.strip())
                                                         for i in line.split('\t')[1:]])
        with open(os.path.join(modelDir, 'function_loss_costs.tsv'), 'r') as f:
            for line in f.readlines()[1:]:
                self.genes[line.split('\t')[0]].functionLossCosts = np.array([float(i.strip())
                                                             for i in line.split('\t')[1:]])

    def test_gene_to_reaction_rules_sensible(self):
        rules = [r.gene_reaction_rule for r in self.model.reactions]
        allRules = ''.join(rules)
        allRules = allRules.replace(' ', '')
        allRules = allRules.replace('and', '')
        allRules = allRules.replace('or', '')
        allRules = allRules.replace('(', '')
        allRules = allRules.replace(')', '')
        for geneName in sorted([g.id for g in self.model.genes], key=len, reverse=True):
            allRules = allRules.replace(geneName, '')
        assert allRules == '', 'gene reaction rules should contain only |gene names|and|or|()|'

    def test_function_loss_equals_gene_loss_in_simple_cases(self):
        assert all([g.old_and_new_costs_identical() for g in self.genes.values() if
                    g.is_simple_single_function()])

    def test_isoenzymes_not_simple_single_function(self):
        assert len([g for g in self.genes.values() if g.is_isoenzyme()
                                        and g.is_simple_single_function()]) == 0





    def test_isoenzyme_pairs_with_only_one_reaction_have_symmetric_costs(self):

        def genes_in_rule(rule):
            """Given a reaction rule, return a list of genes.

            Args:
                rule (str): the reaction rule.

            Returns:
                list(str): the genes.

            """
            genes = set(rule.replace('and', '').replace('or', '').replace('(', '').replace(')', '').split())
            if len(genes) == 0:
                raise UserWarning('ERROR: no genes found in reaction rule.')
            return genes

        isozymesInOneReaction = {gene.name: gene for gene in self.genes.values()
                                 if gene.is_isoenzyme() and gene.number_reactions() == 1}
        isoSimplePairs = set()
        for g in isozymesInOneReaction.values():
            genesInRule = genes_in_rule(g.reactionRules[0])
            if len(genesInRule) == 2:
                if all([i in isozymesInOneReaction for i in genesInRule]):
                    isoSimplePairs.add(tuple(sorted(genesInRule)))
        nIsoSimplePairs = len(isoSimplePairs)
        assert nIsoSimplePairs > 10, 'Should be at least a few simple pairs of isoenzymes.'
        countOldEqual, countNewEqual, countOldZero = 0, 0, 0
        for i, j in isoSimplePairs:
            if np.array_equal(self.genes[i].geneLossCosts, self.genes[j].geneLossCosts):
                countOldEqual += 1
            if np.array_equal(self.genes[i].functionLossCosts, self.genes[j].functionLossCosts):
                countNewEqual += 1
            if np.array_equal(self.genes[i].geneLossCosts, np.zeros(self.genes[i].geneLossCosts.shape)):
                countOldZero += 1
            if np.all(np.isclose(self.genes[j].geneLossCosts, np.zeros(self.genes[j].geneLossCosts.shape), atol = 1e-5)):
                countOldZero += 1
        assert countOldEqual == nIsoSimplePairs
        assert countNewEqual == nIsoSimplePairs
        assert countOldZero == nIsoSimplePairs
