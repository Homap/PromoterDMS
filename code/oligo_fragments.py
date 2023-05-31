import sys
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
# We need two functions for primer design: For PCR of wild-type fragments and another for assembly primer design
def design_primer(fragment, length):
    """
    Design primer for each wildtype fragment
    It slides by one base and produce a new primer pair and assess the best possible pair
    """
    fprimer = fragment[0:length]
    rev_compliment = revcomp(fragment)
    rprimer = rev_compliment[0:length]
    return([fprimer, rprimer])
#*********************************************************************************************
# Plasmid homology arms
HO_locus = plasmid_dict['pGGAselect-YFPins'][412:571]
YFP_locus = plasmid_dict['pGGAselect-YFPins'][658:848]
#*********************************************************************************************
# Function to output mutated and wild-type fragments
def oligo_design(promoter_sequence, fragment_len=150):
    num_oligo_fragments = round(len(promoter_sequence)/fragment_len)
    mutated_oligo_dict = {}
    WT_oligo_dict = {}
    for oligo in range(1, num_oligo_fragments+1):
        if oligo == 1:
            first_fragment = promoter_sequence[0:fragment_len+20]
            key = "oligo_" + str(oligo)
            mutated_oligo_dict[key] = first_fragment
            first_fragment_WT = promoter_sequence[fragment_len:len(promoter_sequence)]
            WT_oligo_dict[key] = first_fragment_WT
        elif oligo == num_oligo_fragments:
            last_fragment = promoter_sequence[((num_oligo_fragments-1)*fragment_len)-20:len(promoter_sequence)]
            key = "oligo_" + str(oligo)
            mutated_oligo_dict[key] = last_fragment
            last_fragment_WT = promoter_sequence[0:((num_oligo_fragments-1)*fragment_len)]
            WT_oligo_dict[key] = last_fragment_WT
        else:
            if num_oligo_fragments > 2:
                mid_fragment = promoter_sequence[((oligo-1)*fragment_len)-20:(oligo*fragment_len)+20]
                key = "oligo_" + str(oligo)
                mutated_oligo_dict[key] = mid_fragment
                mid_fragment_WT_1 = promoter_sequence[0:(oligo-1)*fragment_len]
                mid_fragment_WT_2 = promoter_sequence[oligo*fragment_len:len(promoter_sequence)]
                WT_oligo_dict[key] = [mid_fragment_WT_1, mid_fragment_WT_2]
            else:
                pass
    return([mutated_oligo_dict, WT_oligo_dict])


promoter_f = open("../data/promoters/promoter_sequences.fasta", "r")
promoter_dict = readfasta(promoter_f)

mutated_WT_oligo_dict = {}
for promoter in promoter_dict.keys():
    mutated_WT_oligo_dict[promoter] = oligo_design(promoter_dict[promoter])



fragments_dict = {}
for promoter in promoter_dict.keys():
    promoter_seq = promoter_dict[promoter]
    fragment_list = fragments(promoter_seq)
    for part in fragment_list:
        if promoter in fragments_dict.keys():
            fragments_dict[promoter].append(part)
        else:
            fragments_dict[promoter] = [part]



for promoter in fragments_dict.keys():
    print(promoter)
    for part in fragments_dict[promoter]:
        print(len(part))


PFY1 = "AGGAGACGTTACTTTGTTTATATATATTAGTATGTACAATCGCAAAGAAATGGAGTGATGACATGTTGTAGTATTTAGTATGAGGTTACTGTGTGGGAGGTTTTACCATGATTTTTGGCGAGAACACGCCATGAAATGTCTTTGTACGAAACTCATTACCCGCATTAATATTTTTTTTCTTTTTAAAGCTCAGTTGACCCTTTCTCATTCCCTTCTTAAAACAACTGTGTGATCCTTGAGAAAAGATAAATTACATACACAACATAAACCCAACTACGATCGCAAATT"

round(len(PFY1)/150)

len(PFY1)/2

# Define different experiments: 
# Experiment 1: M1 - W2
# Experiment 2: W1 - M2

# Wild-type PCR primers for PFY1
pfy1_f1 = PFY1[0:144]
pfy1_f2 = PFY1[(144-20):]

PFY1_primers_plus_HO_arm_5 = HO_locus[-30:] + design_primer(pfy1_f1, 20)[0]
PFY1_primers_plus_HO_arm_3 = design_primer(pfy1_f1, 20)[1]
PFY1_primers_plus_YFP_arm_5 = design_primer(pfy1_f2, 20)[0]
PFY1_primers_plus_YFP_arm_3 = (design_primer(pfy1_f2, 20)[1][::-1] + revcomp(YFP_locus, reverse=False)[:30])[::-1]

# Wild-type PCR primers for PFY1
pfy1_f1 = PFY1[0:144-20]
pfy1_f2 = PFY1[(145-20-20):]

PFY1_primers_plus_HO_arm_5 = HO_locus[-30:] + design_primer(pfy1_f1, 20)[0]
PFY1_primers_plus_HO_arm_3 = design_primer(pfy1_f1, 19)[1]
PFY1_primers_plus_YFP_arm_5 = design_primer(pfy1_f2, 20)[0]
PFY1_primers_plus_YFP_arm_3 = (design_primer(pfy1_f2, 20)[1][::-1] + revcomp(YFP_locus, reverse=False)[:30])[::-1]




dsr_PFY1 = Dseqrecord(PFY1)


VMA7 = "TTTCCAATTTTTGATGCGTTCATACCAAGTTTGCTGTTTTTCTTAGTTATTATCATCATCTTTCCCCTTATTGCTGCTGAAACTTGCGTACGTGCGTACTAAAAAAATTTCACTGCGCTTACGCCGGTAAAGGGAAGTAGGAAAACGAAAGAAAAAAAAAAAATAAAGGCGAACGATAACAGGAAAACCAACGTGAATTGCAAGCACTACATTTATTATATCGAGTTGACAAAGAATTAAGAAAATAGAGCACTACATTTCTATATTCATAGCTAGTCTCACTGACGCATAGTAACTAAATC"

round(len(VMA7)/150)
len(VMA7)/2

# Write a function for this

# Wild-type PCR primers for PFY1
vma7_f1 = VMA7[0:(151)]
vma7_f2 = VMA7[(151-20):]

VMA7_primers_plus_HO_arm_5 = HO_locus[-30:] + design_primer(vma7_f1, 21)[0]
VMA7_primers_plus_HO_arm_3 = design_primer(vma7_f1, 20)[1]
VMA7_primers_plus_YFP_arm_5 = design_primer(vma7_f2, 20)[0]
VMA7_primers_plus_YFP_arm_3 = (design_primer(vma7_f2, 20)[1][::-1] + revcomp(YFP_locus, reverse=False)[:30])[::-1]

# Wild-type PCR primers for PFY1
pfy1_f1 = PFY1[0:144-20]
pfy1_f2 = PFY1[(145-20-20):]


OST1 = "GTTAACTTTCTGTGTGTACTTTGAAGTAGCAATTCTTATCGTAATTGTACACCGTCTCGCAACTTGTTTAAAATACTGTTAAAGATTGCATTTCAATGCTTTTCCCTTTTTTGCCACCGGATTCCGTTGCGAGTGGAAAAGTGAAAATGACAAAATTCAAAAGAAATCTATCAATATAACAAAAACGCCAAATCATCAAAAATGATCGGTGTTGAAGAGCTTGGCACTCTTAAAGGCGGCTAGTTCATAATGCTATAGGAATTGACTCAAGGAAACGGACTATGTCTTGTACTGAATACTGTCTTCAATTGCCCATAGAATGCGTTAGTGTTACTCTTCTTCGCGAGCAGGGTAAGTATGCGCGTAATGTTTTTTATTTTCTGAAAGGTTCAAAGTATGCAGAACAAATTAATGTTTGCTTTTTATTAAGCAAGCTACTTCTTTGACAAGTACCCGATTGCTTCTTTAGCTCGGAACAAGACGCAAACTACAAAATATTGGTGCTGAAAAA"

num_fragments = round(len(OST1)/150)
fragment_length = len(OST1)/num_fragments
# Three fragments, two of 170 bps and one of 171 bp
# Three fragments are 171+(170*2)

# Wild-type PCR primers for PFY1
ost1_f1 = OST1[170-30:]
ost1_f2 = OST1[0:170+30]
ost1_f3 = OST1[(170*2)-30:]
ost1_f4 = OST1[:(170*2)+30]

ost1_f1_primer_5 = design_primer(ost1_f1, 20)[0]
ost1_f1_primer_3 = design_primer(ost1_f1, 20)[1] + YFP_locus[:30]

ost1_f2_primer_5 = HO_locus[-30:] + design_primer(ost1_f2, 20)[0]
ost1_f2_primer_3 = design_primer(ost1_f2, 20)[1] 

ost1_f3_primer_5 = design_primer(ost1_f3, 20)[0] 
ost1_f3_primer_3 = design_primer(ost1_f3, 20)[1] + YFP_locus[:30]

ost1_f4_primer_5 = HO_locus[-30:] + design_primer(ost1_f4, 20)[0] 
ost1_f4_primer_3 = design_primer(ost1_f4, 20)[1] 

primers_plus_HO_arm_5 = HO_locus[-30:] + design_primer(vma7_f1, 20)[0]
primers_plus_HO_arm_3 = design_primer(vma7_f1, 20)[1]
primers_plus_YFP_arm_5 = design_primer(vma7_f2, 20)[0]
primers_plus_YFP_arm_3 = design_primer(vma7_f2, 20)[1] + YFP_locus[:30]

#*********************************************************************************************
def mutate_promoter(promoter):
    """
    Function to produce 1 mutation to each of three possibile bases
    on the promoter sequence
    """
    #promoter_dict = readfasta(promoter) # read promoter fasta into a dictionary
    promoter_dict = promoter
    mutated_promoter_dict = {}
    for gene_name in promoter_dict.keys():
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

# Mutated fragment F1
pfy1_f1_mutated = HO_locus[-30:] + PFY1[0:144] # 124 to 144 are bases for the homology arm
# Wild-Type fragment F1
pfy1_f1_WT = PFY1[144-20:] + YFP_locus[:30]

# Mutated fragment F2
pfy1_f2_mutated = PFY1[(145-20-20):]
# Wild-type fragment F1
pfy1_f2_WT = PFY1[0:144-20]

PFY1_promoter = {"PFY1_F1" : "AGGAGACGTTACTTTGTTTATATATATTAGTATGTACAATCGCAAAGAAATGGAGTGATGACATGTTGTAGTATTTAGTATGAGGTTACTGTGTGGGAGGTTTTACCATGATTTTTGGCGAGAACACGCCATGAAATGTCTTTG", 
                 "PFY1_F2" : "CCATGATTTTTGGCGAGAACACGCCATGAAATGTCTTTGTACGAAACTCATTACCCGCATTAATATTTTTTTTCTTTTTAAAGCTCAGTTGACCCTTTCTCATTCCCTTCTTAAAACAACTGTGTGATCCTTGAGAAAAGATAAATTACATACACAACATAAACCCAACTACGATCGCAAATT"}

PFY1_mutated_promoters = mutate_promoter(PFY1_promoter)
pfy1_mutated_oligs = open("mutated_pfy1.txt", "w")
for sequence in PFY1_mutated_promoters.keys():
    if sequence.split("_")[1] == "F1":
        pfy1_mutated_oligs.write(sequence+"\t"+HO_locus[-30:]+PFY1_mutated_promoters[sequence]+"\n")
    elif sequence.split("_")[1] == "F2":
        pfy1_mutated_oligs.write(sequence+"\t"+PFY1_mutated_promoters[sequence]+YFP_locus[:30]+"\n")

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
def split_promoter(mutated_promoter_dict, recog_dict, fixed_len = True, oligo_len, start, max_len, min_len):
    """
    Function to split the promoter sequence into fragments
    We can split the promoter into fixed fragment length but strict length size could not be useful since
    the sequence might not be optimal in terms of the future assembly method. It is best to specify a minimum and maximum
    for the length of the fragment with a start position and then check the sequence for compatibility for recognition site.
    Also evaluate   
    """
    # note!: Adjust the script for when one of the pieces is very short.
    mutated_fragment_dict = {}
    if fixed_len:
        for mutation in mutated_promoter_dict.keys():
            if not mutation in recog_dict.keys():
                fragments = chunks(mutated_promoter_dict[mutation], oligo_len)
                for index, item in enumerate(fragments):
                    # Think about few more bases
                    #fragment_dict = {"fragment1" : RNR1[0:237], "fragment2" : RNR1[237-4:2*237], "fragment3" : RNR1[(2*237)-4:3*237], "fragment4" : RNR1[(3*237)-4:4*237]}
                    fragment_ID = mutation + "_F" + str(index+1)
                    mutated_fragment_dict[fragment_ID] = item
    else:
        min_end = start+min_len
        max_end = start+max_len
        for end in range(min_end, max_end):
            new_oligo_len = 
        oligo_range = [[start, min_end], [start, max_end]]

    return(mutated_fragment_dict)
#*********************************************************************************************
# def oligo_fragments(mutated_fragment_dict):

#*********************************************************************************************
# Split promoter into mutated fragment and wildtype fragment
#*********************************************************************************************

# Gibson assembly is a protocol for joining DNA fragments in-vitro by tratment with a mixture of T5 exconuclease,
# DNA polymerase and Taq DNA ligase. Gibson assembly requires 20-40 bp of perfect homology between 3' and 5' ends for 
# fragments to be joined. The T5 exonuclease chews back each fragment in the 5'-3' direction so that the remaining 3' 
# single stranded overhangs can anneal. Gaps are filled and nicks sealed by polymerase and ligase. 

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