"""
Various utilities used by EnsembleTR
"""

import math

def GetEHScore(conf_invs, CNs, ru_len):
    r"""
    Compute a confidence score for EH genotypes
    
    Parameters
     ----------
    conf_invs : str
        FORMAT/REPCI field from EH
    CNs : str
        FORMAT/REPCN field from EH
    ru_len : int
        Repeat unit length (bp)

    Returns
    -------
    score : float
        Confidence score. 0=low, 1=high
    """
    conf_invs = conf_invs.split("/")
    CNs = CNs.split("/")
    try:
        CNs = [(int(CN) * ru_len) for CN in CNs]
    except:
        print(f"Invalid copy numbers or intervals: {CNs}, {conf_invs}")
        return 0
    score1 = CalcEHAlleleScore(conf_invs[0], CNs[0])
    score2 = CalcEHAlleleScore(conf_invs[1], CNs[1])
    return 0.8 * min(score1, score2) + 0.2 * max(score1, score2)

def CalcEHAlleleScore(conf_inv, allele):
    r"""
    Compute an allele-specific score for
    EH genotypes

    Parameters
    ----------
    conf_inv : str
       low-high (based on one allele in FORMAT/REPCI field)
    allele : int
       Inferred allele (based on one allele in FORMAT/REPCN field)

    Returns
    -------
    score : float
       Confidence score for the allele. 0=low, 1=high
    """
    conf_inv = conf_inv.split("-")
    dist = abs(int(conf_inv[0]) - int(conf_inv[1]))
    if dist > 100:
         return 0
    if allele == 0:
         return 1/math.exp(4 * (dist))
    return 1/math.exp(4 * (dist) / int(allele))
