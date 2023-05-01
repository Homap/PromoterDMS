from Bio.Restriction import *
from Bio.Seq import Seq
from revcomp import revcomp
from fasta import readfasta
from chunks import chunks

"""
This script describes design of oligonucleotides to be synthesized by 
TwistBiosciences with the goal of assembling the fragments using the
Golden Gate Assembly method.
"""

"""
We need to ideally also simulate the PCR of sequences differing by 1 base 
To see the challenges that exist in doing that.
"""

"""
Start with the simplest situation but in the final function, the best solution
is to get to determine the size of the oligonucleotide as well.
"""

#*******************************s
# Nucleotide dictionary
#*******************************
nuc_dict = {'A': ['C', 'G', 'T'],
            'T': ['C', 'G', 'A'],
            'C': ['G', 'A', 'T'],
            'G': ['C', 'A', 'T']}

#*********************************************************************************************
# Specify the inputs
#*********************************************************************************************
promoter = open(sys.argv[1], "r")
plasmid = open(sys.argv[2], "r")

#*********************************************************************************************
# Definition of functions 
#*********************************************************************************************
def mutate_promoter(promoter):
    """
    Function to produce 1 mutation to each of three possibile bases
    on the promoter sequence
    """
    promoter_dict = readfasta(promoter) # read promoter fasta into a dictionary
    gene_name = list(promoter_dict.keys())[0]
    mutated_promoter_dict = {}
    mutated_promoter_dict[gene_name+"_0"+"_WT"] = promoter_dict[gene_name] # initate the mutated dictionary with the wild type sequence
    #for position in range(0, len(promoter_dict[gene_name])):
    for index, item in enumerate(list(promoter_dict[gene_name])):
        nuc = promoter_dict[gene_name][index]
        if nuc == 'A':
            for base in nuc_dict['A']:
                mutated_seq_key = gene_name + "_" + str(index+1) + "_" + nuc + base
                promoter_base_list = list(promoter_dict[gene_name])
                promoter_base_list[index] = base
                mutated_seq = "".join(promoter_base_list)
                mutated_promoter_dict[mutated_seq_key] = mutated_seq
        if nuc == 'T':
            for base in nuc_dict['T']:
                mutated_seq_key = gene_name + "_" + str(index+1) + "_" + nuc + base
                promoter_base_list = list(promoter_dict[gene_name])
                promoter_base_list[index] = base
                mutated_seq = "".join(promoter_base_list)
                mutated_promoter_dict[mutated_seq_key] = mutated_seq
        if nuc == 'C':
            for base in nuc_dict['C']:
                mutated_seq_key = gene_name + "_" + str(index+1) + "_" + nuc + base
                promoter_base_list = list(promoter_dict[gene_name])
                promoter_base_list[index] = base
                mutated_seq = "".join(promoter_base_list)
                mutated_promoter_dict[mutated_seq_key] = mutated_seq
        if nuc == 'G':
            for base in nuc_dict['G']:
                mutated_seq_key = gene_name + "_" + str(index+1) + "_" + nuc + base
                promoter_base_list = list(promoter_dict[gene_name])
                promoter_base_list[index] = base
                mutated_seq = "".join(promoter_base_list)
                mutated_promoter_dict[mutated_seq_key] = mutated_seq
    return(mutated_promoter_dict)
#*********************************************************************************************
def scan_recog_motif(mutated_promoter_dict, RE):
    """
    Scans sequence to see if it contains restriction enzyme
    recognition motif.
    """
    recog_dict = {}
    for mutation in mutated_promoter_dict.keys():
        sequence = Seq(mutated_promoter_dict[mutation])
        recog_site = RE.search(sequence)
        if len(recog_site) > 0:
            if not mutation in recog_dict.keys():
                recog_dict[mutation] = [recog_site]
            else:
                recog_dict[mutation].append(recog_site)
        else:
            print("No recognition site was detected")
    return(recog_dict)
#*********************************************************************************************
def split_promoter(mutated_promoter_dict, recog_dict, oligo_len):
    """
    Function to split the promoter sequence into fragments
    """
    # note!: Adjust the script for when one of the pieces is very short.
    mutated_fragment_dict = {}
    for mutation in mutated_promoter_dict.keys():
        if not mutation in recog_dict.keys():
            fragments = chunks(mutated_promoter_dict[mutation], oligo_len)
            for index, item in enumerate(fragments):
                # Think about few more bases
                #fragment_dict = {"fragment1" : RNR1[0:237], "fragment2" : RNR1[237-4:2*237], "fragment3" : RNR1[(2*237)-4:3*237], "fragment4" : RNR1[(3*237)-4:4*237]}
                fragment_ID = mutation + "_F" + str(index+1)
                mutated_fragment_dict[fragment_ID] = item
    return(mutated_fragment_dict)
#*********************************************************************************************
def oligo_fragments(mutated_fragment_dict):




# The first and last ligo must produce homology with the plasmid sequence after it is digested by BsaI
oligo1 = 20*"N" + "GGTCTC" + "N" + "ACAT" + fragment1 + "N" + "GAGACC" + 20*"N"
oligo2 = 20*"N" + "GGTCTC" + "N" + fragment2 + "N" + "GAGACC" + 20*"N"
oligo3 = 20*"N" + "GGTCTC" + "N" + fragment3 + "N" + "GAGACC" + 20*"N"
oligo4 = 20*"N" + "GGTCTC" + "N" + fragment4 + "TTCA" + "N" + "GAGACC" + 20*"N"

print(mutated_seq_dict)

# Golden Gate assembly with BsaI enzyme


# Bring in plasmid sequence which is a circular DNA





#revcomp(RNR1, reverse=False)

RNR1 = "TCACTGTGTATAGTAGATCTGCAGCTCAGTCACATGAGACATTAACGTATATAAATATATAAATATATATATATCCATTAAATAACGAAGTAGGACATGGTCTAAAAATAGTGCTTCTCATGTATTGATTTTGAGTAGCCATAAGAATAATTATTGATTTCAAGATCTTTGATCTCCTTACATCCAAAATTATAATACCACTTTATGGTATCACGACATCACCACAACAGAGAAGGATTGCCGCATTTGGATATCGTAAACAAAGGCGTTACCATAGAAATGTACTGATTGGCAGAATTACTCTTCAGGAGAATCTTTCATACAAAGGAATTCCATTGGGGAAAATCTCGTTACCAAGTCAATGCTGAACTTTCTATGGCCTTTGTTTACTATCGTTAATTATTTTACGACCACTTCTGGGTAGAAATATTTCGTAGCCCTGGAACGAGCTTGTTTACGCGTTTTATCCCATTATATGGCACCCAAATCAAATTTAAAAAGAAAAAACGCGTAAACAGTGTCGGGTAAGTTCATCCTCTGTTACTTTAATTGCTTCTTTTTTTGAAATTCTAAGTAAACGCGTCATTTTGATCCTCAGGACACAGAAATCCTTGCAGAATCTTATTGGGTGTTGAATAGAGGACGCGTAAAAACGATATGGAAATTTTTTTCATATAGTGTAGAAAGAATAGGTTGGCGTAGGTAGTTTCGTGTTTGATAGAAACCTCCAACAAAGTCTGCAGCTCACGTTTTAGAATAACAAGTTTAGAGTTTATCTTGTTGCCTTTGTTAAGTCAGTACCATTGAATAAAAATTATATAAAGGAGCTAATATTTCATTGTTGGAAAATTACTCTACCATAATTGAAGCATATCTCATCCTTTTCATCCTTTTCAACGCAAGAGAGACACCAACGAACAACACTTTATTTGTTGATATATTAACATC"

fragment1 = RNR1[0:237] 
fragment2 = RNR1[237-4:2*237]
fragment3 = RNR1[(2*237)-4:3*237]
fragment4 = RNR1[(3*237)-4:4*237]

#k = fragment1 + fragment2 + fragment3 + fragment4



oligo_list = [oligo1, oligo2, oligo3, oligo4]

# restriction batch
rb = RestrictionBatch([BsaI])

for oligo in oligo_list:
	Analong = Analysis(rb, Seq(oligo))
	Analong.print_as('map')
	Analong.print_that()

# Now mutate each fragment, create the oligo and check if any of the mutations create additional BsaI recognition motif

# For each nucleotide, mutate it to the other 3

promoter_dict = {"RNR1" : "TCACTGTGTATAGTAGATCTGCAGCTCAGTCACATGAGACATTAACGTATATAAATATATAAATATATATATATCCATTAAATAACGAAGTAGGACATGGTCTAAAAATAGTGCTTCTCATGTATTGATTTTGAGTAGCCATAAGAATAATTATTGATTTCAAGATCTTTGATCTCCTTACATCCAAAATTATAATACCACTTTATGGTATCACGACATCACCACAACAGAGAAGGATTGCCGCATTTGGATATCGTAAACAAAGGCGTTACCATAGAAATGTACTGATTGGCAGAATTACTCTTCAGGAGAATCTTTCATACAAAGGAATTCCATTGGGGAAAATCTCGTTACCAAGTCAATGCTGAACTTTCTATGGCCTTTGTTTACTATCGTTAATTATTTTACGACCACTTCTGGGTAGAAATATTTCGTAGCCCTGGAACGAGCTTGTTTACGCGTTTTATCCCATTATATGGCACCCAAATCAAATTTAAAAAGAAAAAACGCGTAAACAGTGTCGGGTAAGTTCATCCTCTGTTACTTTAATTGCTTCTTTTTTTGAAATTCTAAGTAAACGCGTCATTTTGATCCTCAGGACACAGAAATCCTTGCAGAATCTTATTGGGTGTTGAATAGAGGACGCGTAAAAACGATATGGAAATTTTTTTCATATAGTGTAGAAAGAATAGGTTGGCGTAGGTAGTTTCGTGTTTGATAGAAACCTCCAACAAAGTCTGCAGCTCACGTTTTAGAATAACAAGTTTAGAGTTTATCTTGTTGCCTTTGTTAAGTCAGTACCATTGAATAAAAATTATATAAAGGAGCTAATATTTCATTGTTGGAAAATTACTCTACCATAATTGAAGCATATCTCATCCTTTTCATCCTTTTCAACGCAAGAGAGACACCAACGAACAACACTTTATTTGTTGATATATTAACATC"}





#        mutated_fragments[mutated_seq_key] = 

# Random mutations, combinations of mutations, a selection of mutations from another file like VCF



# restriction batch
rb = RestrictionBatch([BsaI])

# Does the mutation create a new mutation site?
k = []
for oligo in mutated_seq_dict.keys():
	k.append(Analong.print_that(None, title='sequence = multi_site\n', s1='\n no site:\n'))
	Analong.print_as('map')
	Analong.print_that()

# Adding the recognition motifs and primers
mutated_oligo = {}
for mutated_fragment in mutated_seq_dict.keys():
    # Must match the plasmid BsaI digestion
    # first_oligo_fragment = 20*"N" + "GGTCTC" + "N" + "ACAT" + fragment1 + "N" + "GAGACC" + 20*"N"
    oligo = 20*"N" + "GGTCTC" + "N" + mutated_seq_dict[mutated_fragment] + "N" + "GAGACC" + 20*"N"
    mutated_oligo[mutated_fragment] = oligo
    # Must match the plasmid BsaI digestion
    # last_oligo_fragment = 20*"N" + "GGTCTC" + "N" + fragment4 + "TTCA" + "N" + "GAGACC" + 20*"N"

for oligo in mutated_oligo.keys():
	Analong = Analysis(rb, Seq(mutated_oligo[oligo]))
	Analong.print_as('map')
	Analong.print_that()
        
from pydna.dseq import Dseq

seq = Dseq("GGATCCAAA","TTTGGATCC",ovhg=0)

from Bio.Restriction import BsaI
a,b = seq.cut(BsaI)

# Does the mutation create ab exrta digestion site?
# To check that, digest the fragment with the restriction site and count the resulting number of fragments. 
# seq = Dseq(mutated_oligo[oligo],revcomp(mutated_oligo[oligo]),ovhg=0)
dr_len = []
for oligo in mutated_oligo.keys():
    seq = Seq(mutated_oligo[oligo])
    digestion_result = len(list(BsaI.catalyse(seq)))
    if digestion_result == 3:
        dr_len.append(digestion_result)
    elif digestion_result > 3:
        print(oligo, BsaI.catalyse(seq))

overhang5 = "".join(list(g[1]))[0:4]
overhang3 = revcomp("".join(list(g[1]))[-4:], reverse = True)

# Assmble the sequences
from pydna.dseq import Dseq

seq = Dseq("GGATCCAAA","TTTGGATCC",ovhg=0)


Seq("NNNNNNNNNNNNNNNNNNNNGGTCTCNCCCTGTCGT").overhang5(dct=None)
    #seq.cut(BsaI)
    #print(a,c)
    

# We also need the plasmid sequence

# oligo1 = 20*"N" + "GGTCTC" + "N" + "ACAT" + fragment1 + "N" + "GAGACC" + 20*"N"
# oligo2 = 20*"N" + "GGTCTC" + "N" + fragment2 + "N" + "GAGACC" + 20*"N"
# oligo3 = 20*"N" + "GGTCTC" + "N" + fragment3 + "N" + "GAGACC" + 20*"N"
# oligo4 = 20*"N" + "GGTCTC" + "N" + fragment4 + "TTCA" + "N" + "GAGACC" + 20*"N"

# Digest the sequence with restriction site


# Check if the overlapping sequences between fragments are unique so that fragments are assembled in correct
# order based on homology


# Check for BsaI recognition site


# mutated_seq_dict = {}
# for key in promoter_dict.keys():
#     for nuc in promoter_dict[key]:
#         if nuc == 'A':
#             for base in nuc_dict['A']:
#                 mutated_seq = promoter_dict[key].replace('A', base)
#                 if not key in mutated_seq_dict:
#                     mutated_seq_dict[key] = [mutated_seq]
#                 else:
#                     mutated_seq_dict[key].append(mutated_seq)
#         if nuc == 'T':
#             for base in nuc_dict['T']:
#                 mutated_seq = promoter_dict[key].replace('T', base)
#                 if not key in mutated_seq_dict:
#                     mutated_seq_dict[key] = [mutated_seq]
#                 else:
#                     mutated_seq_dict[key].append(mutated_seq)
#         if nuc == 'C':
#             for base in nuc_dict['C']:
#                 mutated_seq = promoter_dict[key].replace('C', base)
#                 if not key in mutated_seq_dict:
#                     mutated_seq_dict[key] = [mutated_seq]
#                 else:
#                     mutated_seq_dict[key].append(mutated_seq)
#         if nuc == 'G':
#             for base in nuc_dict['G']:
#                 mutated_seq = promoter_dict[key].replace('G', base)
#                 if not key in mutated_seq_dict:
#                     mutated_seq_dict[key] = [mutated_seq]
#                 else:
#                     mutated_seq_dict[key].append(mutated_seq)

#print(RNR1)
#print(len(RNR1)*"|")
#print(revcomp(RNR1, reverse=False))

# vcftools --gzvcf 1011Matrix.gvcf.gz --freq --chr chromosome10 --from-bp 434649 --to-bp 435163 --out OST1
# vcftools --gzvcf 1011Matrix.gvcf.gz --freq --chr chromosome7 --from-bp 527327 --to-bp 527632 --out VMA7
# vcftools --gzvcf 1011Matrix.gvcf.gz --freq --chr chromosome15 --from-bp 552882 --to-bp 553175 --out PFY1
# vcftools --gzvcf 1011Matrix.gvcf.gz --freq --chr chromosome10 --from-bp 338004 --to-bp 338270 --out TDH1
# vcftools --gzvcf 1011Matrix.gvcf.gz --freq --chr chromosome4 --from-bp 410056 --to-bp 411828 --out GDP1
# vcftools --gzvcf 1011Matrix.gvcf.gz --freq --chr chromosome12 --from-bp 439824 --to-bp 440470 --out STM1
# vcftools --gzvcf 1011Matrix.gvcf.gz --freq --chr chromosome10 --from-bp 391647 --to-bp 392406 --out RNR2
# vcftools --gzvcf 1011Matrix.gvcf.gz --freq --chr chromosome10 --from-bp 454679 --to-bp 456239 --out TDH2
# vcftools --gzvcf 1011Matrix.gvcf.gz --freq --chr chromosome5 --from-bp 298002 --to-bp 298953 --out RNR1
# vcftools --gzvcf 1011Matrix.gvcf.gz --freq --chr chromosome7 --from-bp 883811 --to-bp 884509 --out TDH3

BsaI_1 = "GGTCTC"
BsaI_2 = "GAGACC"

# Primers cannot contain the restriction site because the primer pairs 
# will then be complementary to each other. 

# First primer: GGTCTCATCACGTATATAGC
# Second primer: CGAACTAATAGACTGAGACC

# 6 bp of recognition motif + 20 bp = 26 bp
# 6 bp of recognition motif + 20 bp = 26 bp
# 300 bp - 52 bp = 248 bp
# 290 bp - 52 bp = 238 bp

# Where to cut? It's not so much of a question, you just need to start
# from one end, cut the sequencing for 248 base pairs and continue to do so.