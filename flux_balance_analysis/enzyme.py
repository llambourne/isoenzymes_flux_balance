import numpy as np


class enzyme:

    tolerance = 4.e-5

    def __init__(self, name):
        self.name = name
        # reaction rules in flux balance model that contain this enzyme
        self.reactionRules = []
        # rank in evolutionary rate
        self.dndsRank = 0
        self.functionLossCosts = []
        self.geneLossCosts = []
        self.blocked = False

    def number_reactions(self):
        if len(self.reactionRules) == 0:
            raise UserWarning('ERROR no reaction rules associated with enzyme')
        return len(self.reactionRules)

    def is_isoenzyme(self):
        """
        Isoenzyme defined as any gene conected to a reaction with
        a rule that contains an 'or'.
        """
        if len(self.reactionRules) == 0:
            raise UserWarning('ERROR no reaction rules associated with enzyme')
        for rule in self.reactionRules:
            if 'or' in rule:
                return True
        return False

    def is_multifunctional(self):
        """Is an enzyme involved in more than one interaction."""
        if len(self.reactionRules) == 0:
            raise UserWarning('ERROR no reaction rules associated with enzyme')
        return len(self.reactionRules) > 1

    def is_simple_single_function(self):
        """
        Is an enzyme involved in only one reaction and that reaction
        only involves that one enzyme.
        """
        if len(self.reactionRules) == 0:
            raise UserWarning('ERROR no reaction rules associated with enzyme')
        return len(self.reactionRules) == 1 and self.reactionRules[0] == self.name

    def old_and_new_costs_identical(self):
        """Are function-loss cost and gene-loss cost equal?

        Use a tolerance for this based on accuaracy we expect from the solver.
        """
        return np.allclose(self.functionLossCosts, self.geneLossCosts, rtol=0., atol=self.tolerance)

    def num_environments_disagree(self):
        """Number of media that have a different value between function-loss
           and gene-loss costs."""
        return sum(np.abs(self.geneLossCosts - self.functionLossCosts)  > self.tolerance)

    def hybrid_cost_1(self):
        """Function loss for multifunctional, gene loss for single-function.

        Pink distribtution in Fig S2.
        """
        if self.is_multifunctional():
            return self.functionLossCosts
        else:
            return self.geneLossCosts

    def hybrid_cost_2(self):
        """Function loss for isoenzymes, gene loss for non-isoenzymes.

        Purple distribution in Fig 2.
        """
        if self.is_isoenzyme():
            return self.functionLossCosts
        else:
            return self.geneLossCosts


    def __str__(self):
        return (self.name + '\n' + 'involved in ' +
                str(self.number_reactions()) +
                ' reactions\nRules of those reactions:\n' +
                str(self.reactionRules))
