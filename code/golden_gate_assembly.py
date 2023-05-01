#!/usr/bin/python
import sys
import re

import fragment_dict as fd

gg_name = open(sys.argv[1], "r")

genotype_dict = {}
for line in gg_name:
    line = line.rstrip().split()
    genotype_dict[line[0]] = line[1:]

print(genotype_dict)

fragment_dict = {}
for genotype in genotype_dict.keys():
    for fragment in genotype_dict[genotype]:
        if not genotype in fragment_dict.keys():
            fragment_dict[genotype] = [fd.gg_fragment_dict[fragment]]
        else:
            fragment_dict[genotype].append(fd.gg_fragment_dict[fragment])

print(fragment_dict)

# AncCerPtdh3-Promoter-Frag1.gb AncCPPtdh3-Promoter-Frag2.gb AncCerPtdh3-Promoter-Frag3.gb AncCPPtdh3-Promoter-Frag4.gb AncCPPtdh3-Promoter-Frag5.gb AncCPPtdh3-Promoter-Frag6.gb GG_assemblies/G9

# fragment_list = [fragment1_seq, fragment2_seq, fragment3_seq, fragment4_seq, fragment5_seq, fragment6_seq]

dg_fragment_dict = {}
for genotype in fragment_dict.keys():
    for fragment in fragment_dict[genotype]:
        motif1 = re.finditer(r"GGTCTC", fragment)
        motif2 = re.finditer(r"GAGACC", fragment)
        if motif1:
            match1 = [match.start() for match in motif1][0]
            # print(match1)
        if motif2:
            match2 = [match.start() for match in motif2][0]
            # print(match2)

        digested_fragment = fragment[match1+7:match2-5]
        if not genotype in dg_fragment_dict.keys():
            dg_fragment_dict[genotype] = [digested_fragment]
        else:
            dg_fragment_dict[genotype].append(digested_fragment)

for genotype in dg_fragment_dict.keys():
    print('>'+genotype)
    print("".join(dg_fragment_dict[genotype]))

# def chunks(s, n):
# 	for start in range(0, len(s), n):
# 		yield s[start:start+n]

# # gg_fasta = open(gg_name, "w")
# # gg_fasta.write(">"+gg_name+"\n")
# # for line in chunks(assembled_fragment, 60):
# # 	gg_fasta.write(line+"\n")




