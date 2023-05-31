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
            first_fragment = promoter_sequence[0:fragment_len+20] # HO_locus[-30:]
            key = "oligo_" + str(oligo)
            mutated_oligo_dict[key] = first_fragment
            first_fragment_WT = promoter_sequence[fragment_len:len(promoter_sequence)]
            WT_oligo_dict[key] = first_fragment_WT
        elif oligo == num_oligo_fragments:
            last_fragment = promoter_sequence[((num_oligo_fragments-1)*fragment_len)-20:len(promoter_sequence)] # YFP_locus[:30]
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

# Gibson assembly is a protocol for joining DNA fragments in-vitro by tratment with a mixture of T5 exconuclease,
# DNA polymerase and Taq DNA ligase. Gibson assembly requires 20-40 bp of perfect homology between 3' and 5' ends for 
# fragments to be joined. The T5 exonuclease chews back each fragment in the 5'-3' direction so that the remaining 3' 
# single stranded overhangs can anneal. Gaps are filled and nicks sealed by polymerase and ligase. 



