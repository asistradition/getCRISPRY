from __future__ import print_function

from lib.bowtie2 import bowtie_make_index
from lib.grna_scorer import find_grna

import argparse
import tempfile
import shutil

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

    args = ap.parse_args()
    with open(args.query, mode="rU") as query_fh:
        queries = list(SeqIO.parse(query_fh, format=args.query_type, alphabet=generic_dna))
        grna_builder(queries, args.genome, genome_type=args.genome_type, outfile_path=args.out, out_num=args.num,
                     cores=args.cpu)


def grna_builder(queries, genome_path, genome_type="fasta", outfile_path=None, out_num=1, cores=1):
    if outfile_path is not None:
        out_fh = open(outfile_path, mode="w")

    try:
        if genome_type != "bowtie":
            bowtie_dir = tempfile.mkdtemp()
            bowtie_idx = bowtie_make_index(genome_path, bowtie_dir)
        else:
            bowtie_idx = genome_path
        data = {}
        for sr in queries:
            candy = []
            try:
                targets = find_grna(sr, bowtie_idx, cores=cores)
            except ValueError:
                continue

            for i in range(out_num):
                try:
                    c_score = int(targets[i][1] * 100) / 100
                    if outfile_path is not None:
                        print("\t".join(map(str, [sr.id, targets[i][0], c_score, targets[i][2]])), file=out_fh)
                    candy.append((targets[i][0], c_score, targets[i][2]))
                except IndexError:
                    candy.append((None, None, None))
            data[sr.id] = candy

            print("{}: Identified {} candidate target sequences".format(sr.id, len(candy)))

        return data
    except:
        raise
    finally:
        if genome_type != "bowtie":
            shutil.rmtree(bowtie_dir)


if __name__ == '__main__':
    grna_build_args()
