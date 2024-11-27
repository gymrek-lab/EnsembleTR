"""
Various utilities used by EnsembleTR
"""

import math
import trtools.utils.utils as utils

IUPAC_map_dict = {'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'],
                  'W': ['A', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['G', 'C', 'T'],
                  'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
                  'N': ['A', 'C', 'G', 'T']}

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


def MotifEquality(motif_1, motif_2):
    r"""
    Compute if two motifs are equal to each other, considering IUPAC characters
    This code computes the list of all possible sequences that each motif can be
    then returns True if there is any overlap between the set of sequences for motif_1
    and motif_2.

    Parameters
    ----------
    motif_1 : str
       first motif
    motif_2 : str
       second motif

    Returns
    -------
    Equality : bool, mutual motif
    """
    if len(motif_1) != len(motif_2):
        return False,None
    potential_sequences_1 = [utils.GetCanonicalMotif(motif) for motif in GetAllPossibleSequences(motif_1)]
    potential_sequences_2 = [utils.GetCanonicalMotif(motif) for motif in GetAllPossibleSequences(motif_2)]
    for seq1 in potential_sequences_1:
        if seq1 in potential_sequences_2:
            return True, seq1
    return False, None


def GetAllPossibleSequences(motif):
    r"""
    Computing all the possible sequences that a motif can be considering IUPAC characters.
    For example, a motif with sequence RGG can be both AGG and GGG. Divide and conquer method is
    used to form all possible sequences.

    Parameters
    ----------
    motif : str
       motif

    Returns
    -------
    All possible sequences : list of str
    """
    possible_seqs = []
    if len(motif) < 1:
        return []
    elif len(motif) == 1:
        if motif[0] in IUPAC_map_dict:
            return IUPAC_map_dict[motif[0]]
        else:
            return motif[0]
    else:
        first_part = GetAllPossibleSequences(motif[0:int(len(motif)/2)])
        second_part = GetAllPossibleSequences(motif[int(len(motif)/2):])
        for seq1 in first_part:
            for seq2 in second_part:
                possible_seqs.append(seq1 + seq2)
    return possible_seqs

