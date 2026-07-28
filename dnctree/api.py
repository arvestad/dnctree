import logging

import fastphylo
from modelmatcher import RateMatrix

import dnctree.exceptions as dncexc
from dnctree.msa import MSA
from dnctree.partialdistancematrix import PartialDistanceMatrix
from dnctree.distances import ml_distance_estimate
from dnctree.simple_algorithm import dnc_tree_simple, dnc_tree_simple_centroid
from dnctree.dnctree_k import dnc_tree_k


dna_models = ["hamming", "jc", "k2p", "tn93"]  # fastphylo's DNA/RNA distance models.

# Protein models fastphylo implements natively (ported from modelmatcher, but
# not a complete copy - e.g. mtArt/MtMAM/MtZoa/FLU and several other
# mitochondrial models aren't there). Anything outside this set falls back to
# dnctree's own modelmatcher + scipy ML estimation.
FASTPHYLO_PROTEIN_MODELS = {
    "BLOSUM62", "DCMUT", "Dayhoff", "HIVb", "HIVw", "JTT", "JTT_DCMut",
    "LG", "MtREV", "PMB", "RtREV", "VT", "WAG", "cpREV",
}


def divide_n_conquer_tree(
    seq_data,
    model_name=None,
    max_n_attempts=100,
    max_clade_size=0.5,
    base_case_size=100,
    first_triple=None,
    algorithm='core-tree'
):
    """
    Input: An object with sequence data for the input alignment,
           and the name of a evolutionary model.
    Returns: a tree as a Newick-formatted string.

    Parameters:
    max_n_attempts: The number of times we try to find a suitable taxa triple to recurse on.
    max_clade_size: The fraction of the taxa that are required before deciding we are happy
                    with a triple. We require that the largest subproblem (clade) has at most this fraction of taxa.
    base_case_size: For clades with at most this many taxa, we just run NJ.
    algorithm:      'core-tree' (default), 'simple' (randomly chooses three taxa to guide setup
                    of three subproblems, joined on a centroid vertex with no distances of its
                    own - see dnc_tree_simple_centroid), or 'simple_old' (the same idea, but
                    joined on a vertex whose distances are estimated - see dnc_tree_simple).
    """
    assert base_case_size >= 3  # This is truly a minimum
    assert max_clade_size >= 0.5  # Pointless to try smaller values
    assert max_n_attempts >= 1  # Truly a minimum

    if isinstance(seq_data, MSA):
        distance_fcn, chosen_model_name = choose_distance_function(
            seq_data.type, model_name
        )
        data_description = f"{seq_data.type.upper()} alignment"
    else:
        datatype = type(seq_data)
        raise ValueError(f'Cannot handle input of type {datatype}')

    dm = PartialDistanceMatrix(seq_data, distance_fcn)  # We will gradually add distances to this "matrix"

    if algorithm == 'simple':
        t = dnc_tree_simple_centroid(dm,
                            seq_data.taxa(),
                            max_n_attempts, # Can be removed?
                            max_clade_size, # Can be removed?
                            base_case_size,
                            first_triple,
                            )
    elif algorithm == 'simple_old':
        t = dnc_tree_simple(dm,
                            seq_data.taxa(),
                            max_n_attempts, # Can be removed?
                            max_clade_size, # Can be removed?
                            base_case_size,
                            first_triple,
                            )
    else:                       # Use Amy's work
        # Using the base_case_size as the core_size.
        t = dnc_tree_k(
            dm, seq_data.taxa(), base_case_size=base_case_size, core_size=base_case_size
        )

    actual_work, fraction_work, n = dm.estimate_computational_savings()
    comment = dm._computational_savings()
    logging.info(f"{comment}")
    aux_info = {
        "distances-computed": actual_work,
        "fraction-computed-distances": round(fraction_work, 3),
        "n-taxa": n,
        "comment": comment,
        "model-name": chosen_model_name,
        "description": data_description,
    }
    if isinstance(seq_data, MSA):
        aux_info["msa-width"] = seq_data.msa_width

    return t, aux_info


def choose_distance_function(seqtype, model_name=None):
    """
    Pick a model depending on sequence type and model_name.
    If the two args are not in harmony (e.g., 'dna' and 'WAG'), throw an error.

    seqtype:    'aa', 'dna', or 'rna'.
    model_name: 'WAG' or any other model implemented in the module "modelmatcher".
                'kimura' is also an acceptable string.
    """
    if seqtype == "aa":
        chosen_model = "WAG"  # Standard choice
        if model_name:  # Replace the standard choice
            if model_name in dna_models:
                raise Exception(
                    "You have requested a DNA model for amino acid sequences."
                )
            else:
                chosen_model = model_name
        if chosen_model in FASTPHYLO_PROTEIN_MODELS:
            distance_fcn = _fastphylo_sequence_distance(chosen_model, fastphylo.SeqType.Protein)
        else:
            model = RateMatrix.instantiate(chosen_model)
            distance_fcn = lambda N: ml_distance_estimate(model, N)  # XXX Correct call?
        return distance_fcn, chosen_model
    elif seqtype in ("dna", "rna"):
        chosen_model = model_name or "k2p"  # fastphylo's own default
        if chosen_model not in dna_models:
            raise dncexc.UnknownModel(
                f"'{chosen_model}' is not a known DNA/RNA model. Choose one of: {dna_models}."
            )
        fp_seq_type = fastphylo.SeqType.RNA if seqtype == "rna" else fastphylo.SeqType.DNA
        distance_fcn = _fastphylo_sequence_distance(chosen_model, fp_seq_type)
        return distance_fcn, chosen_model
    else:
        raise dncexc.UnknownModel(f"Not a valid sequence type: {seqtype}.")


def _fastphylo_sequence_distance(model, fp_seq_type):
    """
    Build a distance function that delegates to fastphylo, computed from the
    two taxa's actual sequences rather than a precomputed count matrix.
    """
    def compute(rec1, rec2):
        aln = fastphylo.Alignment([
            fastphylo.Sequence("t1", str(rec1.seq), seq_type=fp_seq_type),
            fastphylo.Sequence("t2", str(rec2.seq), seq_type=fp_seq_type),
        ])
        dm = fastphylo.distance_matrix(aln, model=model)
        return dm[0, 1]

    compute.operates_on_sequences = True  # Tells PartialDistanceMatrix to skip count_pairs().
    return compute
