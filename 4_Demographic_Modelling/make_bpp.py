#make bpp input
from Bio import SeqIO
import glob

with open("hares_bpp.txt", "w") as out:
    for fasta_file in sorted(glob.glob("hares_*.fasta")):
        seqs = list(SeqIO.parse(fasta_file, "fasta"))
        n_indv = len(seqs)             # number of samples
        length = len(seqs[0].seq)      # sequence length
        out.write(f"{n_indv} {length}\n\n")   # writes in "4 20" format
        for record in seqs:
            out.write(f"^{record.id}\t{str(record.seq)}\n")
        out.write("\n")
