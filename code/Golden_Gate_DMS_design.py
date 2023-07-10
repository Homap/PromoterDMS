from Bio.Seq import Seq
from Bio.Restriction import BsaI
from fasta import readfasta

#******************************************************************
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
#******************************************************************


# fasta_f = open("../data/promoters/promoter_sequences.fasta")
# fasta_dict = readfasta(fasta_f)
# scan_recog_motif(fasta_dict, BsaI)

def goldengate_oligo_design(promoter_sequence, fragment_len=150):
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
