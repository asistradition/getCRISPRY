[![Build Status](https://travis-ci.org/asistradition/getCRISPRY.svg?branch=master)](https://travis-ci.org/asistradition/getCRISPRY)

# getCRISPRY 

This is a python package which facilitates CRISPR guide RNA design. 
It takes sequences of arbitrary length (like a gene) and searches for potential CRISPR guides within the input.
The potential CRISPR guides are then scored for on-target efficency, and for off-target efficiency against a provided
reference genome.

Currently, the only implemented CRISPR enzyme is the spCas9 from pyrogenes. The off-target scoring model is largely from
Doench et al 2016 (https://www.ncbi.nlm.nih.gov/pubmed/26780180)

# Usage

Potential guide sequences are identified by exhaustive search, and then each potential guide is screened for some basic
characteristics, and aligned to the reference genome with bowtie2. The alignments are used to score for off-target activity.
The guide sequences are then ranked and outputted to a TSV file.

# Requirements

- Python 3.X
- Biopython 1.67
- Bowtie2 2.2.6
