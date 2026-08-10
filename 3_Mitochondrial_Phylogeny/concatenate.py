import os
import glob
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Fixed gene order
gene_order = [
    "cob", "nad1", "nad2", "cox1", "cox2",
    "atp8", "atp6", "cox3", "nad3", "nad4l",
    "nad4", "nad5", "nad6"
]

rrn_order = ["rrnS", "rrnL"]

# Regex to detect duplicates like gene_0, gene_1, gene_2
dup_pattern = re.compile(r"(?P<gene>[a-z0-9]+)_(\d+)$", re.IGNORECASE)

# Log file
with open("duplicate_log.txt", "w") as log_file:

    for fasta_file in glob.glob("*_annotated.fas"):
        sample_name = os.path.basename(fasta_file).replace("_annotated.fas", "")
        seq_dict = {}
        rrn_dict = {}
        duplicates = {}

        for record in SeqIO.parse(fasta_file, "fasta"):
            # Split the description by semicolons and whitespace to catch the gene name at the end
            parts = record.description.strip().split(";")
            gene_field = parts[-1].strip().lower()  # e.g., "atp6_0"

            # Check protein genes
            for gene in gene_order:
                if gene in gene_field:
                    match = dup_pattern.match(gene_field)
                    if match and match.group("gene").lower() == gene:
                        # It's a duplicate
                        duplicates.setdefault(gene, []).append(gene_field)
                        # Keep only _0 for concatenation
                        if match.group(2) == "0" and gene not in seq_dict:
                            seq_dict[gene] = str(record.seq)
                    else:
                        # No suffix, just the main gene
                        if gene not in seq_dict:
                            seq_dict[gene] = str(record.seq)

            # Check rRNA genes
            for rrn in rrn_order:
                if rrn.lower() in gene_field:
                    match = dup_pattern.match(gene_field)
                    if match and match.group("gene").lower() == rrn.lower():
                        duplicates.setdefault(rrn, []).append(gene_field)
                        if match.group(2) == "0" and rrn not in rrn_dict:
                            rrn_dict[rrn] = str(record.seq)
                    else:
                        if rrn not in rrn_dict:
                            rrn_dict[rrn] = str(record.seq)

        # Concatenate protein genes
        concatenated_protein = "".join([seq_dict[gene] for gene in gene_order if gene in seq_dict])
        if concatenated_protein:
            protein_record = SeqRecord(Seq(concatenated_protein),
                                       id=sample_name,
                                       description="")
            SeqIO.write(protein_record, f"{sample_name}_allgenes.fas", "fasta")

        # Concatenate rRNAs
        concatenated_rrn = "".join([rrn_dict[rrn] for rrn in rrn_order if rrn in rrn_dict])
        if concatenated_rrn:
            rrn_record = SeqRecord(Seq(concatenated_rrn),
                                   id=sample_name,
                                   description="")
            SeqIO.write(rrn_record, f"{sample_name}_rrn.fas", "fasta")

        # Write duplicates to log
        for gene, dup_list in duplicates.items():
            # Sort duplicates numerically if they have _number
            def sort_key(x):
                m = dup_pattern.match(x)
                return int(m.group(2)) if m else -1
            sorted_dups = sorted(dup_list, key=sort_key)
            log_file.write(f"[{sample_name}] Duplicates for {gene}: {', '.join(sorted_dups)}\n")

