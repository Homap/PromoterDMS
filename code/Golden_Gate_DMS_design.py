from Bio.Seq import Seq
from Bio.Restriction import BsaI
from fasta import readfasta

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