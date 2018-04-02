from __future__ import print_function

from lib.bowtie2 import bowtie_make_index
from lib.grna_scorer import GetCRISPRs

import argparse
import tempfile
import shutil
import multiprocessing

from Bio import SeqIO
from Bio.Alphabet import generic_dna


def grna_build_args():
    ap = argparse.ArgumentParser(description="Search genes for CRISPR target sequences")

    ap.add_argument("--query", dest="query", help="FILE with query genes", metavar="FILE", required=True)
    ap.add_argument("--query_type", dest="query_type", help="Query sequence file TYPE", metavar="TYPE", default="fasta")
    ap.add_argument("--genome", dest="genome", help="Genome sequence file / database", metavar="FILE", required=True)
    ap.add_argument("--genome_type", dest="genome_type", help="Genome sequence file TYPE (fasta, bowtie)",
                    metavar="TYPE", default="fasta")
    ap.add_argument("--out", dest="out", help="Output FILE location", metavar="FILE", default="output.tsv")
    ap.add_argument("--num", dest="num", help="NUMBER of target sequences to save", metavar="NUMBER", default=1,
                    type=int)
    ap.add_argument("--cpu", dest="cpu", help="NUMBER of CPU cores to use", metavar="NUMBER", default=1,
                    type=int)
    ap.add_argument("--location", dest="loc", help="FILE containing the query gene locations", metavar="FILE",
                    default=None)


    args = ap.parse_args()
    with open(args.query, mode="rU") as query_fh:
        queries = list(SeqIO.parse(query_fh, format=args.query_type, alphabet=generic_dna))
        grna_builder(queries, args.genome, genome_type=args.genome_type, loc_path=args.loc, outfile_path=args.out,
                     out_num=args.num, cores=args.cpu)


def grna_builder(queries, genome_path, genome_type="fasta", loc_path= None, outfile_path=None, out_num=1, cores=1):
    if outfile_path is not None:
        out_fh = open(outfile_path, mode="w")

    if loc_path is not None:
        locations = _parse_locations(loc_path)
    else:
        locations = None

    try:
        if genome_type != "bowtie":
            bowtie_dir = tempfile.mkdtemp()
            bowtie_idx = bowtie_make_index(genome_path, bowtie_dir)
        else:
            bowtie_idx = genome_path
        data = {}

        finder = GetCRISPRs(bowtie_idx, locations=locations)
        proc_pool = multiprocessing.Pool(processes=cores, maxtasksperchild=50)

        for sr, candy in proc_pool.imap_unordered(finder.find_grna, queries):
            data[sr.id] = candy

            for i in range(out_num):
                try:
                    c_score = int(candy[i][1] * 100) / 100
                    if outfile_path is not None:
                        print("\t".join(map(str, [sr.id, candy[i][0], c_score, candy[i][2]])), file=out_fh)
                except IndexError:
                    pass

            print("{}: Identified {} candidate target sequences".format(sr.id, min(len(candy), out_num)))

        return data

    except:
        raise

    finally:
        if genome_type != "bowtie":
            try:
                shutil.rmtree(bowtie_dir)
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
