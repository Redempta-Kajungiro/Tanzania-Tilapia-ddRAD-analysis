#!/bin/bash
# populations_stacks.sh
# Parameters:

# --min-samples-per-pop : minimum proportion of individuals in a population required to retain a locus
# --min-maf : minimum minor allele frequency to keep SNP
# -t : number of threads

populations \
  -P ./stacks_output \
  -M ./populations_tilapia.txt \
  -t 8 \
  --min-samples-per-pop 0.75 \
  --min-maf 0.05 \
  --write-single-snp \
  --vcf --structure --fstats \
  -O ./populations_output

