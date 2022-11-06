#!/usr/bin/env python3

import sys
import argparse
from Bio import Seq, SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline


def align_fastas(multifasta_file, alignment_file):
    clustalw_cline = ClustalwCommandline("clustalw", infile=multifasta_file, outfile=alignment_file,
                                         align=True, output="FASTA")
    clustalw_cline()


def translate_snp(dna_file, index, snp):
    coding_dna = SeqIO.read(dna_file, "fasta")

    print("Introducing a SNP, nucleotide {} on position {} will become a {}\n".format(coding_dna[index], index, snp))

    as_list = list(coding_dna.seq)
    as_list[index] = snp

    seq_obj = Seq.Seq("".join(as_list))

    return Seq.translate(seq_obj).strip("*")


def process_alignment(args):
    alignment_file = "/".join(args.multifasta.name.split("/")[:-1]) + "/alignment.fasta"

    align_fastas(args.multifasta.name, alignment_file)

    alignment = AlignIO.read(alignment_file, "fasta")
    mutated_protein = translate_snp(args.codingseq.name, args.index, args.snp)

    # AA count for each position
    conservation_dicts = list()

    for i in range(alignment.get_alignment_length()):
        temp = dict()
        unique_aas = set(alignment[:,i])
        for aa in unique_aas:
            temp[aa] = alignment[:,i].count(aa) / len(alignment[:,i])
        conservation_dicts.append(temp)

    # Scoring
    protein_index = (args.index // 3) + 1

    mutated_aa = mutated_protein[protein_index]
    scoring_dict = conservation_dicts[protein_index]
    
    print("The mutated amino acid {} on postion {} will be scored using this conservation table: {}".format(mutated_aa, protein_index, scoring_dict))

    if mutated_aa in scoring_dict.keys():
        if scoring_dict[mutated_aa] == max(scoring_dict.values()):
            score = 1 - max(scoring_dict.values())
        else:
            score = max(scoring_dict.values()) - scoring_dict[mutated_aa]
    else:
        score = max(scoring_dict.values())

    print("The SNP scored a deleterious value of {}, based on conservation in the protein family. Score range from 0 to 1 (lower is better).".format(score))


def main(args):
    process_alignment(args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate the severity scores for SNPs")

    parser.add_argument("-f", "--multifasta", type=argparse.FileType("r"), help="Path to multiprotein fasta (remove the header)", required=True)
    parser.add_argument("-s", "--codingseq", type=argparse.FileType("r"), help="Path to coding dna fasta", required=True)
    parser.add_argument("-i", "--index", type=int, help="Location for the snp", required=True)
    parser.add_argument("-p", "--snp", type=str, help="Polymorhism on index location", required=True)

    args = parser.parse_args()

    sys.exit(main(args))
