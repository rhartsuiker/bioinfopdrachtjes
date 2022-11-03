#!/usr/bin/env python3

import sys
from Bio import Seq, SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline

def align_fastas(amino_fasta, out_file):
    clustalw_cline = ClustalwCommandline("clustalw", infile=amino_fasta, outfile=out_file,
                                         align=True, output="FASTA")
    clustalw_cline()


def process_alignment(alignment_fasta, nuc_fasta):
    alignment = AlignIO.read(alignment_fasta, "fasta")

    print(alignment)

    # AA count for each position
    scoring_dicts = list()

    for i in range(alignment.get_alignment_length()):
        temp = dict()
        unique_aas = set(alignment[:,i])
        for aa in unique_aas:
            temp[aa] = alignment[:,i].count(aa) / len(alignment[:,i])
        temp = {k: v for k, v in sorted(temp.items(), key=lambda item: item[1], reverse = True)}
        scoring_dicts.append(temp)

    # Scoring
    mutated = nuc_to_amino(nuc_fasta).strip("*")
    
    cost = 0
    cost_per_pos = list()
    
    for aa, d in zip(mutated, scoring_dicts):
        keys = list(d.keys())
        
        pre_cost = cost
        
        if aa != keys[0]:
            cost += d[keys[0]]
            if len(d) > 1:
                for key in keys:
                    if aa == key:
                        cost -= d[key]
                        
        cost_per_pos.append(cost - pre_cost)

    for i,e in enumerate(cost_per_pos):
        print("Position {} has a cost of {}".format(i,e))
    print("The total cost for the mutated sequence is: {}".format(cost))


def nuc_to_amino(nuc_fasta):
    coding_dna = SeqIO.read(nuc_fasta, "fasta")
    return Seq.translate(coding_dna.seq)


def main():
    multi_fasta = "/homes/rhartsuiker/bioinf3/hbb/hbb_protein_fam.fasta"
    alignment_file = "/homes/rhartsuiker/bioinf3/hbb/hbb_protein_aln.fasta"
    nuc_fasta = "/homes/rhartsuiker/bioinf3/hbb/coding_seq.fasta"

    align_fastas(multi_fasta, alignment_file)

    process_alignment(alignment_file, nuc_fasta)
    
# ~ def main(args):
    # ~ parser = argparse.ArgumentParser(description="Calculate the severity scores for SNPs"
                                                 # ~ "in a MSA.")

    # ~ parser.add_argument("file", metavar="F", type=str,
                        # ~ help="The file path to the Multiple sequence alignment, only multi Fasta files")

    # ~ parser.add_argument("-l", metavar="location", type=str,
                        # ~ help="A single location for the SNP to calculate the severity")


    # ~ parser.add_argument("-s", metavar="save", type=str, choices=["amino", "alignment", "all"],
                        # ~ help="What the program wil save and write to the destination denoted "
                             # ~ "by -out, you can choose from the following options:"
                             # ~ "amino; only saves the amino translation,"
                             # ~ "alignment; only saves the alignment information,"
                             # ~ "all; saves amino & alignment information."
                             # ~ "If no option was given, will not save data and continue to only "
                             # ~ "calculate the single SNP location.")

    # ~ parser.add_argument("-out", metavar="output", type=str,
                        # ~ help="The path and filename for saving the desired information, "
                             # ~ "make sure to denote it as an .csv file.")


if __name__ == "__main__":
    sys.exit(main())
