from __future__ import print_function

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

DNA_EQUALITY = {"A": ["A", "N", "R", "M", "W", "V", "D", "H"],
                "T": ["T", "N", "Y", "K", "W", "B", "D", "H"],
                "G": ["G", "N", "R", "K", "K", "B", "D", "V"],
                "C": ["C", "N", "Y", "M", "K", "B", "V", "H"],
                "N": ["A", "T", "G", "C", "N", "R", "Y", "M", "K", "W", "V", "D", "H", "B"],
                "R": ["R", "A", "G"],
                "Y": ["Y", "C", "T"],
                "M": ["M", "A", "C"],
                "S": ["G", "C", "S"],
                "W": ["W", "A", "T"],
                "K": ["G", "T", "K"],
                "V": ["V", "A", "G", "C"],
                "D": ["D", "A", "G", "T"],
                "H": ["H", "A", "T", "C"],
                "B": ["B", "T", "G", "C"]}


def str_compare_degenerate(str1, str2, alphabet=DNA_EQUALITY):
    if len(str1) != len(str2):
        return False

    for i in range(len(str1)):
        if str1[i].upper() in alphabet[str2[i].upper()]:
            continue
        else:
            return False

    return True


def fix_seq_object_list(obj_list, alphabet=generic_dna, make_type="str"):
    new_list = []

    for obj in obj_list:
        new_list.append(fix_seq_object(obj, alphabet=alphabet, make_type=make_type))

    return new_list


def fix_seq_object(obj, alphabet=generic_dna, make_type="str"):
    if isinstance(obj, Seq):
        if make_type == "seq":
            return obj
        obj = SeqRecord(obj, id="", description="")
    elif isinstance(obj, str):
        if make_type == "str":
            return obj
        obj = SeqRecord(Seq(obj, alphabet=alphabet), id="", description="")
    elif isinstance(obj, SeqRecord):
        if make_type == "seqrecord":
            return obj
        pass
    else:
        raise ValueError("Not a string, Bio.Seq, or Bio.SeqRecord object")

    if make_type == "seqrecord":
        return obj
    elif make_type == "seq":
        return obj.seq
    elif make_type == "str":
        return str(obj.seq)
    else:
        raise AttributeError("make_type must be str, seq, or seqrecord")

def rc_string(seq, alphabet=generic_dna):
    seq = fix_seq_object(seq, make_type="str")
    return str(Seq(seq, alphabet=generic_dna).reverse_complement())
