#!/usr/bin/env python3

import sys
from Bio import Seq, SeqIO, AlignIO
from Bio.Align.Applications import ClustalwCommandline

def align_fastas(amino_fasta, out_file):
    clustalw_cline = ClustalwCommandline("clustalw", infile=amino_fasta, outfile=out_file,
                                         align=True, output="FASTA")
    clustalw_cline()

def process_alignment(alignment_fasta, ):
    alignment = AlignIO.read(alignment_fasta, "fasta")
    
    print(alignment)

    # AA count for each position
    pos_aa_count = list()

    for i in range(alignment.get_alignment_length()):
        temp = dict()
        unique_aas = set(alignment[:,i])
        for aa in unique_aas:
            temp[aa] = alignment[:,i].count(aa)
        pos_aa_count.append(temp)

    # Percentage conservation in each position
    aa_conserve = [max(d.values()) / len(alignment[:,0]) * 100 for d in pos_aa_count]

    conserv_percent = list()
    
    for i,c in enumerate(aa_conserve):
        print("pos {} has {}% conservation".format(i,int(c)))
        
    # Scoring
    # ~ nuc_to_amino()
    
def nuc_to_amino(nuc_fasta):
    coding_dna = SeqIO.read(nuc_fasta, "fasta")
    return Seq.translate(coding_dna.seq)


def main():
    nuc_to_amino("/homes/rhartsuiker/bioinf3/hbb/coding_seq.fasta")

    fasta_file = "/homes/rhartsuiker/bioinf3/hbb/hbb_protein_fam.fasta"
    alignment_file = "/homes/rhartsuiker/bioinf3/hbb/hbb_protein_aln.fasta"

    align_fastas(fasta_file, alignment_file)

    process_alignment(alignment_file)

if __name__ == "__main__":
    sys.exit(main())
    
# ~ AlignIO.write(alignment, "/homes/rhartsuiker/bioinf3/hbb/hbb_protein_fam.phy", "phylip")
# ~ AlignIO.write(alignment, "/homes/rhartsuiker/bioinf3/hbb/hbb_protein_fam.sth", "stockholm")
