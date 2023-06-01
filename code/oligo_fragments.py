import sys
from Bio.Restriction import *
from Bio.Seq import Seq
from revcomp import revcomp
from fasta import readfasta
from chunks import chunks

"""
This script describes design of oligonucleotides with the goal of 
assembling the fragments using the Gibson Assembly method.
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
        if oligo == 1: # Each oligo contains overlap and the region to be mutated
            #first_fragment = promoter_sequence[0:fragment_len+20] # HO_locus[-30:]
            first_fragment = promoter_sequence[0:fragment_len]
            first_fragment_overlap = promoter_sequence[fragment_len:fragment_len+20]
            key_m = "oligo_" + str(oligo) + "_M"
            key_w = "oligo_" + str(oligo) + "_W"
            mutated_oligo_dict[key_m] = {"HO_arm": HO_locus[-30:], "oligo": first_fragment, "overlap": first_fragment_overlap}
            first_fragment_WT = promoter_sequence[fragment_len:len(promoter_sequence)]
            WT_oligo_dict[key_w] = {"oligo": first_fragment_WT, "YFP_arm": YFP_locus[:30]}
        elif oligo == num_oligo_fragments:
            #last_fragment = promoter_sequence[((num_oligo_fragments-1)*fragment_len)-20:len(promoter_sequence)] # YFP_locus[:30]
            last_fragment = promoter_sequence[((num_oligo_fragments-1)*fragment_len):len(promoter_sequence)]
            last_fragment_overlap = promoter_sequence[((num_oligo_fragments-1)*fragment_len)-20:((num_oligo_fragments-1)*fragment_len)]
            key_m = "oligo_" + str(oligo) + "_M"
            key_w = "oligo_" + str(oligo) + "_W"
            mutated_oligo_dict[key_m] = {"overlap": last_fragment_overlap, "oligo": last_fragment, "YFP_arm": YFP_locus[:30]}
            last_fragment_WT = promoter_sequence[0:((num_oligo_fragments-1)*fragment_len)]
            WT_oligo_dict[key_w] = {"HO_arm": HO_locus[-30:], "oligo": last_fragment_WT}
        else:
            if num_oligo_fragments > 2:
                mid_fragment = promoter_sequence[((oligo-1)*fragment_len)-20:(oligo*fragment_len)+20]
                mid_fragment_overlap1 = promoter_sequence[((oligo-1)*fragment_len)-20:((oligo-1)*fragment_len)]
                mid_fragment_overlap2 = promoter_sequence[(oligo*fragment_len):(oligo*fragment_len)+20]
                key_m = "oligo_" + str(oligo) + "_M"
                key_w = "oligo_" + str(oligo) + "_W"
                mutated_oligo_dict[key_m] = {"overlap1": mid_fragment_overlap1, "oligo": mid_fragment, "overlap2": mid_fragment_overlap2}
                mid_fragment_WT_1 = promoter_sequence[0:(oligo-1)*fragment_len]
                mid_fragment_WT_2 = promoter_sequence[oligo*fragment_len:len(promoter_sequence)]
                WT_oligo_dict[key_w] = {"HO_arm": HO_locus[-30:], "oligo_1": mid_fragment_WT_1,"oligo_2": mid_fragment_WT_2, "YFP_arm": YFP_locus[:30]}
            else:
                pass
    return([{"mutated": mutated_oligo_dict}, {"wildtype": WT_oligo_dict}])
#*********************************************************************************************
def mutate_oligo(oligo_seq):
    """
    Function to produce 1 mutation to each of three possibile bases
    on the promoter sequence
    """
    mutated_promoter_list = []
    for index, item in enumerate(list(oligo_seq)):
        nuc = list(oligo_seq)[index]
        if nuc == 'A':
            for base in nuc_dict['A']:
                promoter_base_list = list(oligo_seq)
                promoter_base_list[index] = base
                mutated_seq = "".join(promoter_base_list)
                mutated_promoter_list.append(mutated_seq)
        elif nuc == 'T':
            for base in nuc_dict['T']:
                promoter_base_list = list(oligo_seq)
                promoter_base_list[index] = base
                mutated_seq = "".join(promoter_base_list)
                mutated_promoter_list.append(mutated_seq)
        elif nuc == 'C':
            for base in nuc_dict['C']:
                promoter_base_list = list(oligo_seq)
                promoter_base_list[index] = base
                mutated_seq = "".join(promoter_base_list)
                mutated_promoter_list.append(mutated_seq)
        elif nuc == 'G':
            for base in nuc_dict['G']:
                promoter_base_list = list(oligo_seq)
                promoter_base_list[index] = base
                mutated_seq = "".join(promoter_base_list)
                mutated_promoter_list.append(mutated_seq)
        else:
            print("Error: Illigal character")
    return(mutated_promoter_list)
#*********************************************************************************************

# Gibson assembly is a protocol for joining DNA fragments in-vitro by tratment with a mixture of T5 exconuclease,
# DNA polymerase and Taq DNA ligase. Gibson assembly requires 20-40 bp of perfect homology between 3' and 5' ends for 
# fragments to be joined. The T5 exonuclease chews back each fragment in the 5'-3' direction so that the remaining 3' 
# single stranded overhangs can anneal. Gaps are filled and nicks sealed by polymerase and ligase. 

promoter_f = open("../data/promoters/promoter_sequences.fasta", "r")
promoter_dict = readfasta(promoter_f)

mutated_WT_oligo_dict = {}
for promoter in promoter_dict.keys():
    mutated_WT_oligo_dict[promoter] = oligo_design(promoter_dict[promoter])


for gene in mutated_WT_oligo_dict.keys():
    for oligo_m in mutated_WT_oligo_dict[gene][0]['mutated'].keys():
        oligo_to_mutate = mutated_WT_oligo_dict[gene][0]['mutated'][oligo_m]['oligo']
        mutated_oligo_list = mutate_oligo(oligo_to_mutate)
        for mutation in mutated_oligo_list:
            print(gene, oligo_m, mutated_WT_oligo_dict[gene][0]['mutated'][oligo_m], mutation)