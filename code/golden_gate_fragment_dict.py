#!/usr/bin/python
import sys

fragment1 = open(sys.argv[1])
fragment2 = open(sys.argv[2])
fragment3 = open(sys.argv[3])
fragment4 = open(sys.argv[4])
fragment5 = open(sys.argv[5])
fragment6 = open(sys.argv[6])
fragment7 = open(sys.argv[7])
fragment8 = open(sys.argv[8])
fragment9 = open(sys.argv[9])
fragment10 = open(sys.argv[10])
fragment11 = open(sys.argv[11])
fragment12 = open(sys.argv[12])

gg_fragment_dict = {}

def gb_to_fasta(fragment):
    sequence = []
    for line in fragment:
        line = line.strip("\n").split()
        if line[0].startswith("LOCUS"):
            header = line[1]
        if line[0].isdigit():
            sequence.append("".join(line[1:]).upper())
    return("".join(sequence))

fragment1_seq = gb_to_fasta(fragment1)
fragment2_seq = gb_to_fasta(fragment2)
fragment3_seq = gb_to_fasta(fragment3)
fragment4_seq = gb_to_fasta(fragment4)
fragment5_seq = gb_to_fasta(fragment5)
fragment6_seq = gb_to_fasta(fragment6)
fragment7_seq = gb_to_fasta(fragment7)
fragment8_seq = gb_to_fasta(fragment8)
fragment9_seq = gb_to_fasta(fragment9)
fragment10_seq = gb_to_fasta(fragment10)
fragment11_seq = gb_to_fasta(fragment11)
fragment12_seq = gb_to_fasta(fragment12)

gg_fragment_dict = {'cer1':fragment1_seq,
                    'cer2':fragment2_seq,
                    'cer3':fragment3_seq,
                    'cer4':fragment4_seq,
                    'cer5':fragment5_seq,
                    'cer6':fragment6_seq,
                    'cp1':fragment7_seq,
                    'cp2':fragment8_seq,
                    'cp3':fragment9_seq,
                    'cp4':fragment10_seq,
                    'cp5':fragment11_seq,
                    'cp6':fragment12_seq
}

print("gg_fragment_dict = "+str(gg_fragment_dict))