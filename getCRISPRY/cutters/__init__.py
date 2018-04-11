from __future__ import print_function

from ..lib.degenerate_tools import str_compare_degenerate, fix_seq_object, rc_string

class Cutter:
    """
    Cutter is an class that holds the framework for a CRISPR-type RNA-directed DNA binding enzyme.
    The class parameters should be set within a child that extends this class.
    Scoring algorithms can also be overloaded within the child class to give cutter-specific scoring
    """

    # These class parameters must be set when Cutter is extended
    cutter_name = None
    PAM = None
    PAM_3prime = True
    spacer_length = None

    # These class parameters can be set when Cutter is extended.
    seq = None
    name = None
    default_penalty = 0.85

    @classmethod
    def __str__(cls):
        if cls.name is not None:
            return cls.name
        else:
            raise AttributeError("Undefined name")

    @classmethod
    def on_target_score(cls, guide_sequence, target_sequence):

        guide_sequence, target_sequence, pam_seq = cls._process_sequences(guide_sequence, target_sequence)

        if str_compare_degenerate(guide_sequence, target_sequence) and str_compare_degenerate(pam_seq, cls.PAM):
            return 1
        else:
            return 0

    @classmethod
    def off_target_score(cls, guide_sequence, target_sequence):
        """
        Calculates the score of a guide sequence and a target sequence using the class parameters for the
        endonuclease

        :param guide_sequence: str
            String corresponding to the guide RNA sequence (without the PAM)
        :param target_sequence: str
            String corresponding to the target site in the genome

        :return: score [numeric]
            Cutting frequency score of the target sequence with the specified guide sequence
        """

        guide_sequence, target_sequence, pam_seq = cls._process_sequences(guide_sequence, target_sequence)
        guide_sequence = str(guide_sequence).replace("T", "U")

        if str_compare_degenerate(pam_seq, cls.PAM):
            score = 1
        else:
            score = 0

        # Walk through the guide sequence and compare it to the target sequence
        # If there is a mismatch, look up the mismatch score
        g_s = list(guide_sequence)
        t_s = list(target_sequence)

        for i, g_chr in enumerate(g_s):
            try:
                t_chr = t_s[i]
            except IndexError:
                t_chr = "-"

            if g_chr is t_chr:
                continue
            elif g_chr is "U" and t_chr is "T":
                continue

            if cls.PAM_3prime:
                score *= cls.default_penalty ** i
            else:
                score *= cls.default_penalty ** (cls.spacer_length - i)

        return score

    @classmethod
    def remove_pam(cls, seq):
        """
        Removes the PAM sequence from a given target. Raises a PAMNotFoundError if the PAM sequence is not present.

        :param seq:
            Sequence to process

        :return seq:
            Sequence with PAM sequence removed
        """
        target_length = cls.spacer_length + len(cls.PAM)
        if len(seq) != target_length:
            raise ValueError(
                "Target sequence length is {}bp [{} takes {}bp]".format(len(seq), cls.cutter_name, target_length))

        if cls.PAM_3prime:
            if not str_compare_degenerate(seq[cls.spacer_length:], cls.PAM):
                raise PAMNotFoundError
            return seq[:cls.spacer_length]
        else:
            if not str_compare_degenerate(seq[:len(cls.PAM)], cls.PAM):
                raise PAMNotFoundError
            return seq[len(cls.PAM):]

    @classmethod
    def acceptable_sequence(cls, seq, stem_length=4, loop_length=6, homopolymer_max=3):
        """
        Checks to see if a sequence is an acceptable guide on structure only (does not look for PAM sequence,
        valid_target should be used for that). Raises a ValueError if the sequence has unsuitable regions.

        Currently looks for hairpins and for homopolymer regions using the laziest possible implementation

        :param seq: str
            Sequence to check

        :return: bool
            Returns True
        """

        seq = fix_seq_object(seq, make_type="str")

        # Check for homopolymer
        a_chr = ""
        chr_count = 0
        for i in range(len(seq)):
            if seq[i] == a_chr:
                chr_count += 1
            else:
                a_chr = seq[i]
                chr_count = 0

            if chr_count > homopolymer_max:
                raise StructureUnsuitedError("Homopolymer")

        # Quick and dirty check for hairpins
        for si, sj in cls._hairpin_slider(seq, stem_length=stem_length, loop_length=loop_length):
            if si == rc_string(sj):
                raise StructureUnsuitedError("Hairpin")

        return True

    @classmethod
    def valid_target(cls, seq, start, stop):
        """
        Checks to see if a given sequence matches the PAM sequence of the cutter. Raises a ValueError if the length of
        the provided sequence is not equal to the length of the PAM sequence. Raises a PAMNotFoundError if the sequence
        does not match the PAM sequence.

        :param seq: str
            Candidate PAM sequence
        :param start: int
            Start position of the candidate sequence
        :param stop:
            Stop position of the candidate sequence

        :return m_st, m_en, strand: (int, int, str)
             Returns location coordinates for the start and end of the genomic target, and the strand "+" or "-"
        """

        seq = fix_seq_object(seq, make_type="str")

        if len(seq) != len(cls.PAM):
            raise ValueError("Sequence length does not match PAM length")

        if str_compare_degenerate(cls.PAM, seq):
            if cls.PAM_3prime:
                m_st = start - cls.spacer_length
                m_en = stop
            else:
                m_st = start
                m_en = stop + cls.spacer_length
            return m_st, m_en, "+"

        elif str_compare_degenerate(cls.PAM, rc_string(seq)):
            if cls.PAM_3prime:
                m_st = start
                m_en = stop + cls.spacer_length
            else:
                m_st = start - cls.spacer_length
                m_en = stop
            return m_st, m_en, "-"

        raise PAMNotFoundError("Sequence {} does not match PAM [{}]".format(seq, cls.PAM))

    @classmethod
    def pam_length(cls):
        return len(cls.PAM)

    @classmethod
    def _process_sequences(cls, guide_sequence, target_sequence):

        # Uppercase the sequence strings
        guide_sequence = str(guide_sequence).upper()
        target_sequence = str(target_sequence).upper()

        # Pull off the PAM sequence from the 5' or 3' end
        if cls.PAM_3prime:
            pam_seq = target_sequence[-1 * len(cls.PAM):]
            target_sequence = target_sequence[0:len(target_sequence) - len(pam_seq)]
        else:
            pam_seq = target_sequence[0:len(cls.PAM)]
            target_sequence = target_sequence[len(cls.PAM):]

        return guide_sequence, target_sequence, pam_seq

    @classmethod
    def _hairpin_slider(cls, s, stem_length=4, loop_length=None):

        if len(s) - stem_length <= 0:
            raise StopIteration

        for start_1 in range(len(s) - stem_length):

            stop_1 = start_1 + stem_length
            seq1 = s[start_1:stop_1]

            for k in range(min(loop_length, len(s[stop_1:]))):

                start_2 = k + stop_1
                stop_2 = start_2 + stem_length

                if stop_2 > len(s):
                    break

                seq2 = s[start_2:start_2 + stem_length]

                yield seq1, seq2


class PAMNotFoundError(LookupError):
    pass


class StructureUnsuitedError(ValueError):
    pass
