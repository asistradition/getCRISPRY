from __future__ import print_function

from .degenerate_tools import fix_seq_object, fix_seq_object_list, rc_string
from .bowtie2 import bowtie_get_hits, bowtie_check_index
from ..cutters.cas9 import cas9
from ..cutters import PAMNotFoundError

from Bio.SeqUtils import GC


class GetCRISPRs:

    def __init__(self, bowtie_idx, at_optimal=0.65, pos_optimal=0.5, cutter=cas9, score_cutoff=0.01, weight=None,
                 cores=1, verbose=False, locations=None):
        """
        Sets the unchanging parameters for the CRISPR target search algorithm.
        Implemented as a class to make mapping the find_grna function to a list of sequences easier

        Required Arguments:

        :param bowtie_idx: str
            Path string to the bowtie sequence index

        Keyword Arguments:

        :param at_optimal: float
            Optimal ratio of AT bp over total
        :param pos_optimal: float
            Optimal position of grna relative to the target region (0.5 is middle, 0 is 5'end, 1 is 3'end)
        :param cutter: cutter
            The CRISPR cutter that will be used (required for spacer length and PAM sequence)
        :param score_cutoff: float
            The minimum target score required to consider a gRNA/target pair as being possible to cut
            Off-targets under this score will be ignored
        :param weight: function
            Weighted guide scoring function taking on-target score, off-target score, off-target count,
            0-optimal AT ratio, and 0-optimal normalized position score
        :param cores: int
            Number of CPU cores to use for search
        :param locations: dict[str] = {(str, int, int, str)}
            Dict of tuple of chromosome ID, start, stop, and strand that is keyed to the query sequence ids.
            If this is not set, bowtie2 will be used to obtain the location in the genome of the target region.
        """

        # Check to make sure the bowtie index exists
        try:
            bowtie_check_index(bowtie_idx)
            self.idx = bowtie_idx
        except FileNotFoundError:
            print("{} bowtie index not found".format(bowtie_idx))
            raise

        self.at_optimal = at_optimal
        self.pos_optimal = pos_optimal
        self.cutter = cutter
        self.score_cutoff = score_cutoff
        self.cores = cores
        self.verbose = verbose
        self.locations = locations

        if weight is None:
            self.weight = self.default_weight
        else:
            self.weight = weight

    def find_grna(self, region_seq):
        """
        Searches through a given sequence block to identify optimal target sequences for CRISPR RNA-directed DNA binding

        :param region_seq: str
            Sequence string containing the target region

        :return processed_candidates: list [str, float, int]
            List of tuples consisting of candidate target sequence (including PAM), weighted score (lower is better), and
            total number of predicted off-target cuts. Returns list sorted by weighted score, so that the best candidate
            is at the 0 position.
        """

        # Find the location of the region to search in the genome
        current_location = None

        if self.locations is not None:
            try:
                current_location = self.locations[region_seq.id]
            except KeyError:
                pass

        if current_location is None:
            if self.verbose:
                print("{}\tIdentifying Genomic Location".format(region_seq.id), end="\r")
            current_location = _find_location(fix_seq_object(region_seq, make_type="str"), self.idx, cores=self.cores)

        # Generate a list of possible sequences
        candidate_target_sequences = get_candidates(region_seq, cutter=self.cutter)

        # Align the possible sequencs to the genome for off-target mapping
        if self.verbose:
            print("{}\tLocating {} targets in genome".format(region_seq.id, len(candidate_target_sequences)), end="\r")
        genome_wide_hits = bowtie_get_hits(
            fix_seq_object_list(candidate_target_sequences.keys(), make_type="seqrecord"),
            self.idx, cores=self.cores)

        processed_candidates = []
        for i, candy_seq in enumerate(candidate_target_sequences.keys()):

            if self.verbose:
                print("{}\t Scanning {} of {}".format(region_seq.id, i + 1, len(candidate_target_sequences)), end="\r")

            try:
                self.cutter.acceptable_sequence(candy_seq)
            except ValueError:
                continue

            on_score = self.cutter.on_target_score(self.cutter.remove_pam(candy_seq), candy_seq)
            off_score = 0
            off_count = 0

            candy_start, candy_stop, candy_strand = candidate_target_sequences[candy_seq].pop()

            for gen_seq, gen_chr, gen_pos, gen_strand, num_mismatch in genome_wide_hits[candy_seq]:

                score = self.cutter.off_target_score(self.cutter.remove_pam(candy_seq), gen_seq)

                if score < self.score_cutoff or current_location is None:
                    continue

                for r_chr, r_start, r_end, _ in current_location:
                    if r_chr != gen_chr or gen_pos < r_start or gen_pos > r_end:
                        off_count += 1
                        break

                if off_count > 0 and off_score == 0:
                    off_score = score
                elif off_count > 0:
                    off_score = max(score * (1 + off_score), off_score * (1 + score))
                else:
                    pass

            candy_at = 1 - GC(self.cutter.remove_pam(candy_seq)) / 100 - self.at_optimal
            candy_pos_norm = self.pos_optimal - ((candy_start + candy_stop) / 2) / len(region_seq)
            candy_weighted_score = self.weight(on_score, off_score, off_count, candy_at, candy_pos_norm)
            processed_candidates.append((candy_seq, candy_weighted_score, off_count))

        processed_candidates = sorted(processed_candidates, key=lambda x: x[1])

        if self.verbose:
            print("{}\t Scanned {} Candidates\tRetained {}".format(region_seq.id,
                                                                   len(candidate_target_sequences),
                                                                   len(processed_candidates)))

        return region_seq, processed_candidates

    def default_weight(self, ons, offs, offn, a_t, p_n):
        return 10 * (1 - ons) + offn * ((1 + offs) ** 3) + (abs(a_t) * 3) + p_n ** 2


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
        """
        Generator for a sliding window

        :param s: str
            String to slide through
        :param window: int
            Window size

        :yield: (str, int, int)
            Window substring, start position, stop position
        """
        for k in range(len(s) - window + 1):
            yield (s[k:k + window], k, k + window)

    # Create a sliding window of the PAM motif length
    for i, target in enumerate(_sliding_window(query_seq, cutter.pam_length())):
        try:
            # Check to see if the current window is a valid target for the cutter. Raise PAMNotFoundError if not.
            seq, k_st, k_end = target
            start, stop, strand = cutter.valid_target(seq, k_st, k_end)

            # Keep moving if the PAM is too close to one end of the query sequence to get a full spacer
            if start < 0 or stop < 0:
                continue
            if start > len(query_seq) or stop > len(query_seq):
                continue

            # Get the sequence of the spacer + PAM from the query sequence
            if strand is "+":
                candy = query_seq[start:stop].upper()
            else:
                candy = rc_string(query_seq[start:stop].upper())

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

