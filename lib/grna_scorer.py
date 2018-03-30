from __future__ import print_function

from lib.degenerate_tools import fix_seq_object, fix_seq_object_list, rc_string
from lib.bowtie2 import bowtie_get_hits, bowtie_check_index
from cutters.cas9 import cas9
from cutters import PAMNotFoundError

from Bio.SeqUtils import GC

DEFAULT_WEIGHT = lambda ons, offs, offn, a_t, p_n: 10 * (1 - ons) + offn * ((1 + offs) ** 3) + (abs(a_t) * 3) + p_n ** 2


def find_grna(region_seq, bowtie_idx, at_optimal=0.65, pos_optimal=0.5, cutter=cas9, score_cutoff=0.01, weight=None,
              cores=1, location=None):
    """
    Searches through a given sequence block to identify optimal target sequences for CRISPR RNA-directed DNA binding

    :param region_seq: str
        Sequence string containing the target region
    :param bowtie_idx: str
        Path string to the bowtie sequence index
    :param at_optimal: float
        Optimal ratio of AT bp over total
    :param pos_optimal: float
        Optimal position of grna relative to the target region (0.5 is middle, 0 is 5'end, 1 is 3'end)
    :param cutter: cutter
        The CRISPR cutter that will be used (required for spacer length and PAM sequence)
    :param score_cutoff: float
        The minimum target score required to consider a gRNA/target pair as being possible to cut
        Off-targets under this score will be ignored
    :param weight: lambda
        Weighted guide scoring function taking on-target score, off-target score, off-target count, 0-optimal AT ratio,
        and 0-optimal normalized position score
    :param cores: int
        Number of CPU cores to use for search
    :param location: [(str, int, int, str)]
        List of tuple of chromosome ID, start, stop, and strand that corresponds to the query sequence location(s) in
        the genome.If this is not set, bowtie will be used to try to obtain it.

    :return processed_candidates: list [str, float, int]
        List of tuples consisting of candidate target sequence (including PAM), weighted score (lower is better), and
        total number of predicted off-target cuts. Returns list sorted by weighted score, so that the best candidate
        is at the 0 position.
    """

    # Check to make sure the bowtie index exists
    try:
        bowtie_check_index(bowtie_idx)
    except FileNotFoundError:
        print("{} bowtie index not found".format(bowtie_idx))
        raise

    # Set the grna scoring algorithm to default if not passed in
    if weight is None:
        weight = DEFAULT_WEIGHT

    # Find the location of the region to search in the genome
    if location is None:
        location = _find_location(fix_seq_object(region_seq, make_type="str"), bowtie_idx, cores=cores)

        # Generate a list of possible sequences
    candidate_target_sequences = get_candidates(region_seq, cutter=cutter)

    # Align the possible sequencs to the genome for off-target mapping
    genome_wide_hits = bowtie_get_hits(fix_seq_object_list(candidate_target_sequences.keys(), make_type="seqrecord"),
                                       bowtie_idx, cores=cores)

    processed_candidates = []
    for candy_seq in candidate_target_sequences.keys():

        try:
            cutter.acceptable_sequence(candy_seq)
        except ValueError:
            continue

        on_score = cutter.on_target_score(cutter.remove_pam(candy_seq), candy_seq)
        off_score = 0
        off_count = 0

        candy_start, candy_stop, candy_strand = candidate_target_sequences[candy_seq].pop()

        for gen_seq, gen_chr, gen_pos, gen_strand, num_mismatch in genome_wide_hits[candy_seq]:

            score = cutter.off_target_score(cutter.remove_pam(candy_seq), gen_seq)

            if score < score_cutoff or location is None:
                continue

            for r_chr, r_start, r_end, _ in location:
                if r_chr != gen_chr or gen_pos < r_start or gen_pos > r_end:
                    off_count += 1
                    break

            if off_count > 0 and off_score == 0:
                off_score = score
            elif off_count > 0:
                off_score = max(score * (1 + off_score), off_score * (1 + score))
            else:
                pass

        candy_at = 1 - GC(cutter.remove_pam(candy_seq)) / 100 - at_optimal
        candy_pos_norm = pos_optimal - ((candy_start + candy_stop) / 2) / len(region_seq)
        candy_weighted_score = weight(on_score, off_score, off_count, candy_at, candy_pos_norm)
        processed_candidates.append((candy_seq, candy_weighted_score, off_count))

    processed_candidates = sorted(processed_candidates, key=lambda x: x[1])
    return processed_candidates


def get_candidates(query_seq, cutter=cas9):
    """
    Takes a sequence and identifies all the sequences which can possibly be used as guides for the specified CRISPR
    cutter

    :param query_seq: sequence
        The sequence to search as a string, Bio.Seq or Bio.SeqRecord object
    :param cutter:  cutter
        The CRISPR cutter that will be used (required for spacer length and PAM sequence)

    :return {seq: [(start, stop, strand)]: {str: [(int, int, str)]
        Returns a dict of candidates keyed by a sequence string with a start, a stop, and a strand
        character (+/-) in a list of tuples
    """

    query_seq = fix_seq_object(query_seq, make_type="str")

    candidates = {}

    def _sliding_window(s, window):
        for k in range(len(s) - window + 1):
            yield (s[k:k + window], k, k + window)

    for i, target in enumerate(_sliding_window(query_seq, cutter.pam_length())):
        try:
            seq, k_st, k_end = target
            start, stop, strand = cutter.valid_target(seq, k_st, k_end)
            if strand is "+":
                candy = query_seq[start:stop].upper()
            else:
                candy = rc_string(query_seq[start:stop].upper())
            if start < 0 or stop < 0:
                continue
            if start > len(query_seq) or stop > len(query_seq):
                continue

            try:
                candidates[str(candy)].append((start, stop, strand))
            except KeyError:
                candidates[str(candy)] = [(start, stop, strand)]

        except PAMNotFoundError:
            pass

    return candidates


def _find_location(query, bowtie_idx, require_unique=False, cores=1):
    region = bowtie_get_hits(fix_seq_object(query, make_type="seqrecord"), bowtie_idx, cores=cores)

    if len(region) == 0 and require_unique:
        raise ValueError("Cannot identify sequence location in genome")
    elif len(region) > 1 and require_unique:
        raise ValueError("Region location indeterminate")

    regions = None

    for qs in region:
        for hit in region[qs]:
            _, r_chr, r_pos, r_strand, _ = hit
            r_start = r_pos
            r_end = r_pos + len(query)
            try:
                regions.append((r_chr, r_start, r_end, r_strand))
            except AttributeError:
                regions = [(r_chr, r_start, r_end, r_strand)]

    return regions
