import numpy as np
import cobra


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


def get_reaction_from_name(reactions, name):
    """Return reaction with name.

    Args:
        reactions (list(reactions)): the reactions.
        name (str): name of the desired reaction.

    Returns:
        reaction: the corresponding reaction.

    """
    matches = [r for r in reactions if r.name == name]
    if len(matches) > 1:
        raise UserWarning('ERROR: duplicate reactions in model.')
    elif len(matches) == 0:
        raise UserWarning('WARNING: could not find reaction: ' + name)
    else:
        return matches[0]


def get_exchange_reactions(model):
    """The input reactions of a model.

    Args:
        model (cobra.model): FBA model.

    Returns:
        list(cobra.reaction): exchange reactions.

    """
    S = model.to_array_based_model().S
    indices = list(np.array(range(len(model.reactions)))[((S != 0.).sum(axis=0) == 1).A1])
    return [model.reactions[i] for i in indices]


def find_media_reactions(sourceNames, exchangeReactions):
    """Map nutrients to exchange reactions in FBA model.

    Args:
        sourceNames (list(str)): names of the metabolites.
        exchangeReactions (list(reactions)): the exchange reactions.

    Returns:
        list(reactions): the corresponding reactions.

    """
    reactions = []
    for nutrient in sourceNames:
        matches = [r for r in exchangeReactions if r.name == nutrient+' exchange']
        if len(matches) > 1:
            raise UserWarning('ERROR: duplicate reactions in model.')
        elif len(matches) == 0:
            print 'WARNING: could not find corresponding reaction for:', nutrient
        else:
            reactions.append(matches[0])
    if len(reactions) != len(sourceNames):
        print 'WARNING: did not find all minimal media reactions'
        print 'found', len(reactions), 'out of', len(sourceNames)
    return reactions


def minimal_media_model(modelPath, carbonSource, nitrogenSource):
    """Load FBA model and set bounds to a minimal media.

    Args:
        modelPath (str): path to FBA model xml file.
        carbonSource (str): name of metabolite.
        nitrogenSource (str): name of metabolite.

    Returns:
        cobra.model: FBA model.

    """
    model = cobra.io.read_sbml_model(modelPath)
    exchangeReactions = get_exchange_reactions(model)
    model.reactions.get_by_id([r.id for r in exchangeReactions
                              if r.name == 'D-glucose exchange'][0]).lower_bound = 0.
    model.reactions.get_by_id([r.id for r in exchangeReactions
                              if r.name == 'ammonium exchange'][0]).lower_bound = 0.
    if carbonSource == nitrogenSource:
        get_reaction_from_name(exchangeReactions,
                               carbonSource + ' exchange').lower_bound = -20.
    else:
        get_reaction_from_name(exchangeReactions,
                               carbonSource + ' exchange').lower_bound = -10.
        get_reaction_from_name(exchangeReactions,
                               nitrogenSource + ' exchange').lower_bound = -10.
    return model


def minimal_media(model, carbonSource, nitrogenSource):
    """Set FBA model to minimal media environment.

    Args:
        model (cobra.model): FBA model.
        carbonSource (str): name of metabolite.
        nitrogenSource (str): name of metabolite.

    Returns:
        cobra.model: FBA model with minimal media envirnoment settings.
    """
    standardNutrients = ['H+', 'iron(2+)', 'oxygen', 'phosphate', 'potassium',
                         'sodium', 'sulphate', 'water']
    stndrdRxns = set([i + ' exchange' for i in standardNutrients])
    exchangeReactions = get_exchange_reactions(model)
    for r in exchangeReactions:
        if r.name in stndrdRxns:
            r.lower_bound = -1000
        if r.name == carbonSource + ' exchange' or r.name == nitrogenSource + ' exchange':
            if carbonSource == nitrogenSource:
                r.lower_bound = -20.
            else:
                r.lower_bound = -10.
        else:
            r.lower_bound = 0.
    return model


def load_sd_minus_his(modelPath):
    """FBA model with conditions that are the same as the
       media used in the Constanza 2009 Science paper

    Args:
        modelPath (str): path to FBA model xml file.

    Returns:
        cobra.model: FBA model.

    """
    #TODO: add a bunch of checks here
    model = cobra.io.read_sbml_model(modelPath)
    # Note: different models have different names for their reactions
    sdMinusHis = set([i+' exchange' for i in ['biotin', 'choline',
                                              'myo-inositol', 'uracil',
                                              'L-alanine', 'L-arginine',
                                              'L-asparagine', 'L-aspartate',
                                              'L-cystein', 'glycine',
                                              'L-glutamate', 'L-glutamine',
                                              'L-isoleucine', 'L-leucine',
                                              'L-lysine', 'L-methionine',
                                              'L-phenylalanine', 'L-proline',
                                              'L-serine', 'L-threonine',
                                              'L-tryptophan', 'L-tyrosine',
                                              'L-valine']])
    exchangeReactions = get_exchange_reactions(model)
    for r in exchangeReactions:
        if r.name in sdMinusHis:
            r.lower_bound = -10.
    return model


def load_media(modelPath, mediaName):
    """FBA model with rich media conditions.

    Args:
        modelPath (str): path to FBA model xml file.
        mediaName (str): type of media to load.

    Returns:
        cobra.model: FBA model.

    """
    sdMinusHis = set(['biotin', 'choline', 'myo-inositol', 'uracil',
                      'L-alanine', 'L-arginine', 'L-asparagine', 'L-aspartate',
                      'L-cystein', 'glycine', 'L-glutamate', 'L-glutamine',
                      'L-isoleucine', 'L-leucine', 'L-lysine', 'L-methionine',
                      'L-phenylalanine', 'L-proline', 'L-serine', 'L-threonine',
                      'L-tryptophan', 'L-tyrosine', 'L-valine'])
    ypd = sdMinusHis.union(set(['L-histodine', 'riboflavin', 'thiamine(1+)',
                                'thymidine', 'nicotinate', '4-aminobenzoate',
                                 '(R)-pantothenate', 'pyridoxine']))
    nutrients = {'sd_minus_his': sdMinusHis,
                 'ypd': ypd}
    if mediaName not in nutrients:
        raise UserWarning(mediaName + ' not available. Choose ' +
                          ' or '.join(nutrients.keys()))
    model = cobra.io.read_sbml_model(modelPath)
    # Note: different models have different names for their reactions
    nutrientRxns = set([i+' exchange' for i in nutrients[mediaName]])
    exchangeReactions = get_exchange_reactions(model)
    for r in exchangeReactions:
        if r.name in sdMinusHis:
            r.lower_bound = -10.
    return model
