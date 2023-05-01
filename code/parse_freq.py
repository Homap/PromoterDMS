#!/usr/bin/python
import sys

frq = open(sys.argv[1], "r")

frq_dict = {}
for line in frq:
    if not line.startswith("CHROM"):
        line = line.strip().split()
        if line[2] == "2":
            allele1 = line[4].split(":")[0]
            allele2 = line[5].split(":")[0]
            if len(allele1) == 1 and len(allele2) == 1:
                key = allele1+"_"+allele2
                print(key + "\t" + line[5].split(":")[1])
        
