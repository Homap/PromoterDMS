import sys
from revcomp import revcomp
from fasta import readfasta
from chunks import chunks

"""
This script describes design of wild-type oligonucleotides with the goal of assembling the fragments using the Gibson Assembly method.
To add: Wild-type fragments here have YFP and HO homology arms. To produce these fragments by PCR, we cannot have that, rather
the plasmid primers should contain promoter-specific homology arms to get linearised and then use those arms for the Gibson assembly.
"""

#*********************************************************************************************
# Specify the inputs
#*********************************************************************************************
promoter_f = open("../data/promoters/promoter_sequences.fasta", "r")
plasmid_f = open("../data/pggaselect-yfpins.fasta", "r")

promoter_dict = readfasta(promoter_f)
plasmid_dict = readfasta(plasmid_f)
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
def gibson_oligo_design(promoter_sequence, fragment_len=150):
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
            key_m = "oligo_" + "last" + "_M"
            key_w = "oligo_" + "last" + "_W"
            mutated_oligo_dict[key_m] = {"overlap": last_fragment_overlap, "oligo": last_fragment, "YFP_arm": YFP_locus[:30]}
            last_fragment_WT = promoter_sequence[0:((num_oligo_fragments-1)*fragment_len)]
            WT_oligo_dict[key_w] = {"HO_arm": HO_locus[-30:], "oligo": last_fragment_WT}
        else:
            if num_oligo_fragments > 2:
                mid_fragment = promoter_sequence[((oligo-1)*fragment_len):(oligo*fragment_len)]
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
# Gibson assembly is a protocol for joining DNA fragments in-vitro by tratment with a mixture of T5 exconuclease,
# DNA polymerase and Taq DNA ligase. Gibson assembly requires 20-40 bp of perfect homology between 3' and 5' ends for 
# fragments to be joined. The T5 exonuclease chews back each fragment in the 5'-3' direction so that the remaining 3' 
# single stranded overhangs can anneal. Gaps are filled and nicks sealed by polymerase and ligase. 

mutated_WT_oligo_dict = {}
for promoter in promoter_dict.keys():
    mutated_WT_oligo_dict[promoter] = oligo_design(promoter_dict[promoter])

# Add coordinate of the mutation + add primers to linearise the plasmid + add primers to build wild type fragments + add primers to build mutated fragments (this is just to test the protocol)
#Gene
#Experiment
#Mutated_fragment_overlap_left
#Mutated_fragment
#Mutated_fragment_overlap_right
#Wild_type_fragment_1
#Wild_type_fragment_2
#Linearizing_primer_1
#Linearizing_primer_2
#Wild_type_primer_1
#Wild_type_primer2
#Mutated_fragment_primer_1
#Mutated_fragment_primer2

HO_reverse_primer = design_primer(HO_locus, 22)[1]
YFP_forward_primer = design_primer(YFP_locus, 22)[0]

gene_list = ['PFY1_YOR122C', 'VMA7_YGR020C', 'OST1_YJL002C', 'STM1_YLR150W', 'TDH3_YGR192C', 'RNR2_YJL026W', 'RNR1_YER070W', 'TDH1_YJL052W', 'TDH2_YJR009C', 'GPD1_YDL022W']

header = ["gene", "experiment", "oligo_M", "M_overlap_left", "M_fragment", "M_overlap_right", "WT_f1", "WT_f2", "Plasmid_primer_HO", "Plasmid_primer_YFP", "WT_f1_primer_forward", "WT_f1_primer_reverse", "WT_f2_primer_forward", "WT_f2_primer_reverse"]
print("\t".join(header))
for gene in gene_list:#mutated_WT_oligo_dict.keys():
    experiment_count = 0
    for oligo_m in mutated_WT_oligo_dict[gene][0]['mutated'].keys():
        experiment_count += 1
        experiment = "Experiment_" + str(experiment_count)
        oligo_base = '_'.join(oligo_m.split("_")[0:2])
        oligo_w = oligo_base + "_W"
        if oligo_m.split("_")[1] == "1":
            HO_arm = mutated_WT_oligo_dict[gene][0]['mutated'][oligo_m]['HO_arm']
            overlap = mutated_WT_oligo_dict[gene][0]['mutated'][oligo_m]['overlap']
            oligo_to_mutate = mutated_WT_oligo_dict[gene][0]['mutated'][oligo_m]['oligo'] + overlap
            oligo_to_mutate_complement = revcomp(oligo_to_mutate, reverse=False, complement=True) 
            seqlist_1 = ["5' ", oligo_to_mutate, " 3'"]
            complist_1 = ["3' ", oligo_to_mutate_complement, " 5'"]
            print("****************************************************************")
            print("".join(seqlist_1)[0:30]+"..."+"".join(seqlist_1)[-30::])
            print("".join(complist_1)[0:30]+"..."+"".join(complist_1)[-30::])
            print("****************************************************************")
            #YFP_arm_W = mutated_WT_oligo_dict[gene][1]['wildtype'][oligo_w]['YFP_arm']
            wildtype_oligo = mutated_WT_oligo_dict[gene][1]['wildtype'][oligo_w]['oligo'] 
            wildtype_oligo_complement = revcomp(wildtype_oligo, reverse=False, complement=True)
            seqlist_2 = ["5' ", wildtype_oligo, " 3'"]
            complist_2 = ["3' ", wildtype_oligo_complement, " 5'"]
            print("****************************************************************")
            print("".join(seqlist_1)[0:30]+"..."+"".join(seqlist_1)[-30::], "".join(seqlist_2)[0:30]+"..."+"".join(seqlist_2)[-30::])
            print("".join(complist_1)[0:30]+"..."+"".join(complist_1)[-30::], "".join(complist_2)[0:30]+"..."+"".join(complist_2)[-30::])
            print("****************************************************************")
            wildtype_oligo_forward = design_primer(wildtype_oligo, 22)[0]
            wildtype_oligo_reverse = design_primer(wildtype_oligo, 22)[1]
            YFP_linearizing_primer = wildtype_oligo[-22:] + design_primer(YFP_locus, 22)[0]
            output_list = [gene, experiment, oligo_m, HO_arm, oligo_to_mutate, overlap, wildtype_oligo, "NA", HO_reverse_primer, YFP_linearizing_primer, wildtype_oligo_forward, wildtype_oligo_reverse, "NA", "NA"]
            #print("\t".join(output_list))
        # elif oligo_m.split("_")[1] == "last":
        #     oligo_to_mutate = mutated_WT_oligo_dict[gene][0]['mutated'][oligo_m]['oligo']
        #     YFP_arm = mutated_WT_oligo_dict[gene][0]['mutated'][oligo_m]['YFP_arm']
        #     overlap = mutated_WT_oligo_dict[gene][0]['mutated'][oligo_m]['overlap']
        #     #HO_arm_W = mutated_WT_oligo_dict[gene][1]['wildtype'][oligo_w]['HO_arm']
        #     wildtype_oligo = mutated_WT_oligo_dict[gene][1]['wildtype'][oligo_w]['oligo']
        #     wildtype_oligo_forward = design_primer(wildtype_oligo, 22)[0]
        #     wildtype_oligo_reverse = design_primer(wildtype_oligo, 22)[1]
        #     HO_linearizing_primer = wildtype_oligo[0:22] + design_primer(HO_locus, 22)[1]
        #     # right and left becomes different depending which oligo is the mutated one.
        #     output_list = [gene, experiment, oligo_m, YFP_arm, oligo_to_mutate, overlap, wildtype_oligo, "NA", HO_linearizing_primer, YFP_forward_primer, wildtype_oligo_forward, wildtype_oligo_reverse, "NA", "NA"]
        #     print("\t".join(output_list))
        # else:
        #     oligo_to_mutate = mutated_WT_oligo_dict[gene][0]['mutated'][oligo_m]['oligo']
        #     overlap1 = mutated_WT_oligo_dict[gene][0]['mutated'][oligo_m]['overlap1']
        #     overlap2 = mutated_WT_oligo_dict[gene][0]['mutated'][oligo_m]['overlap2']
        #     HO_arm_W = mutated_WT_oligo_dict[gene][1]['wildtype'][oligo_w]['HO_arm']
        #     YFP_arm_W = mutated_WT_oligo_dict[gene][1]['wildtype'][oligo_w]['YFP_arm']

        #     wildtype_oligo1 = mutated_WT_oligo_dict[gene][1]['wildtype'][oligo_w]['oligo_1']
        #     wildtype_oligo2 = mutated_WT_oligo_dict[gene][1]['wildtype'][oligo_w]['oligo_2']
        #     wildtype_oligo1_forward = design_primer(wildtype_oligo1, 22)[0]
        #     wildtype_oligo1_reverse = design_primer(wildtype_oligo1, 22)[1]
        #     wildtype_oligo2_forward = design_primer(wildtype_oligo2, 22)[0]
        #     wildtype_oligo2_reverse = design_primer(wildtype_oligo2, 22)[1]
        #     HO_linearizing_primer = wildtype_oligo[0:22] + design_primer(HO_locus, 22)[0]
        #     YFP_linearizing_primer = wildtype_oligo[-22:] + design_primer(YFP_locus, 22)[0]
        #     output_list = [gene, experiment, oligo_m, overlap1, oligo_to_mutate, overlap2, wildtype_oligo1, wildtype_oligo2, HO_linearizing_primer, YFP_linearizing_primer, wildtype_oligo1_forward, wildtype_oligo1_reverse, wildtype_oligo2_forward, wildtype_oligo2_reverse]
        #     print("\t".join(output_list))

# To do:
# Check if adding the two fragments make the complete fragment