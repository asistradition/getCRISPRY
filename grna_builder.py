from __future__ import print_function

from lib.bowtie2 import bowtie_make_index
from lib.grna_scorer import GetCRISPRs
from cutters.cas9 import cas9

import argparse
import tempfile
import shutil
import multiprocessing

from Bio import SeqIO
from Bio.Alphabet import generic_dna

CUTTER_MAP = {"cas9": cas9}


def grna_build_args():
    ap = argparse.ArgumentParser(description="Search genes for CRISPR target sequences")

    ap.add_argument("--query", dest="query", help="FILE with query genes", metavar="FILE", required=True)
    ap.add_argument("--query_type", dest="query_type", help="Query sequence file TYPE", metavar="TYPE", default="fasta")
    ap.add_argument("--genome", dest="genome", help="Genome sequence file / database", metavar="FILE", required=True)
    ap.add_argument("--genome_type", dest="genome_type", help="Genome sequence file TYPE (fasta, bowtie)",
                    metavar="TYPE", default="fasta")
    ap.add_argument("--out", dest="out", help="Output FILE location", metavar="FILE", default="output.tsv")
    ap.add_argument("--num", dest="num", help="NUMBER of target sequences to save", metavar="NUMBER", default=None,
                    type=int)
    ap.add_argument("--cpu", dest="cpu", help="NUMBER of CPU cores to use", metavar="NUMBER", default=1,
                    type=int)
    ap.add_argument("--location", dest="loc", help="FILE containing the query gene locations", metavar="FILE",
                    default=None)
    ap.add_argument("--cutter", dest="cutter",
                    help="The NAME of the supported CRISPR-type enzyme to identify sites for",
                    metavar="NAME", default="cas9")

    args = ap.parse_args()

    # Turn the string corresponding to the CRISPR enzyme into a Cutter object
    try:
        cutter = CUTTER_MAP[args.cutter.lower()]
    except KeyError:
        cutter = None
        print("CRISPR Cutting Enzyme {} Not Supported".format(args.cutter))
        print("Supported Enzymes Are: ", end="")
        print("\t".join(CUTTER_MAP.keys()))
        exit(1)

    # Parse the input sequences file into a list of SeqRecords
    with open(args.query, mode="r") as query_fh:
        queries = list(SeqIO.parse(query_fh, format=args.query_type, alphabet=generic_dna))

    # Parse the input locations file into a dict of location lists, keyed by ID
    if args.loc is not None:
        locations = _parse_locations(args.loc)
    else:
        locations = None

    # Pass everything onto the main guide finder function
    grna_builder(queries, args.genome, genome_type=args.genome_type, locations=locations, outfile_path=args.out,
                 out_num=args.num, cores=args.cpu, cutter=cutter)


def grna_builder(queries, genome_path, genome_type="fasta", locations=None, outfile_path=None, out_num=None, cores=1,
                 cutter=cas9):
    """
    Find guide RNA sequences that will work well in a provided genome for a list of SeqRecord sequences

    Required Arguments:

    :param queries: [SeqRecord]
        Sequence regions to identify CRISPR guide RNAs for
    :param genome_path: str
        A path to a genome file - either FASTA or a bowtie2 index

    Keyword Arguments:

    :param genome_type: str
        The type of genome file - either "fasta" or "bowtie"
    :param locations: {[(str, int, int, str)]}
        A dict of lists of locations (chr, start, stop, strand), keyed by sequence ID
    :param outfile_path: str
        A path to an output file to save sequences as a TSV
    :param out_num: int
        The number of sequences per region to save. Keeps all if set to None
    :param cores: int
        The number of cores to use. Turned into the number of subprocesses in a mp.pool
    :param cutter: Cutter
        A CRISPR Cutter object

    :return data: {(str, float, int)}
        A dict, keyed by sequence ID, of lists of guide RNAs, sorted so that the most optimal is first
    """
    if outfile_path is not None:
        out_fh = open(outfile_path, mode="w")

    try:

        # Create a temporary bowtie index if the genome was provided as a fasta file
        if genome_type.lower() != "bowtie":
            bowtie_dir = tempfile.mkdtemp()
            bowtie_idx = bowtie_make_index(genome_path, bowtie_dir)
        else:
            bowtie_idx = genome_path

        # Create a CRISPR search object and a process pool
        finder = GetCRISPRs(bowtie_idx, locations=locations, cutter=cutter)
        proc_pool = multiprocessing.Pool(processes=cores, maxtasksperchild=50)

        # Map the CRISPR find_grna method across the list of seqrecord objects
        data = {}
        for sr, candy in proc_pool.imap_unordered(finder.find_grna, queries):

            # Save the candidate sequences into the data dict keyed by ID
            if out_num is not None:
                data[sr.id] = candy[:min(out_num, len(candy))]
            else:
                data[sr.id] = candy

            print("{}: Identified {} candidate target sequences".format(sr.id, len(data[sr.id])))

            # Print the candidate sequencs into the output TSV if outfile_path is set
            if outfile_path is not None:
                for grna in data[sr.id]:
                    print("\t".join(map(str, [sr.id, grna[0], int(grna[1] * 100) / 100, grna[2]])), file=out_fh)

        return data

    except:
        raise

    finally:
        # Clean up the temp file and close any outstanding file handles
        if genome_type.lower() != "bowtie":
            try:
                shutil.rmtree(bowtie_dir)
            except UnboundLocalError:
                pass
        if outfile_path is not None:
            try:
                out_fh.close()
            except UnboundLocalError:
                pass


def _parse_locations(loc_path):
    locations = {}
    with open(loc_path, mode="rU") as loc_fh:
        for line in loc_fh:
            l_arr = line.strip().split()
            try:
                locations[l_arr[0]].append((l_arr[1], int(l_arr[2]), int(l_arr[3]), l_arr[4]))
            except KeyError:
                try:
                    locations[l_arr[0]] = [(l_arr[1], int(l_arr[2]), int(l_arr[3]), l_arr[4])]
                except IndexError:
                    continue
            except IndexError:
                continue
    return locations


if __name__ == '__main__':
    grna_build_args()
