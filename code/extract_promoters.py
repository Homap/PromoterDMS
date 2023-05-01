#!/usr/bin/python
import sys
from fasta import readfasta
import random

# python3.9 extract_promoters.py ../data/yeast_genome/GCF_000146045.2_R64_genomic.fna gene_coordinates.txt

genome = open(sys.argv[1], "r")
coord = open(sys.argv[2], "r")

# Read genome fasta into a dictionary
genome_dict = readfasta(genome)

# Read promoter coordinates into dictionary
coord_dict = {}
for line in coord:
    if not line.startswith("gene"):
        line = line.strip().split()
        refseq = line[4]
        gene = line[0]
        start = int(line[1])-1
        end = int(line[2])
        key = gene + ":" + refseq
        coord_dict[key] = sorted([start, end])

# Print out promoter sequences into a file
for key in coord_dict:
    gene = key.split(":")[0]
    chrom = key.split(":")[1]
    start = coord_dict[key][0]
    end = coord_dict[key][1]
    promoter = genome_dict[chrom][start:end].upper()
    GC = ((promoter.count("G") + promoter.count("C"))/len(promoter))*100
    print(">"+gene+"_"+str(len(promoter))+"_"+str(len(promoter)/248))
    print(promoter)


