#!/usr/bin/env python3
"""
Generate the LSD2 date file (dates.txt) for the IQ-TREE dating analysis.

Usage:
    python3 make_dates_file.py hares_mt_allgenes_final_aligned_withAncient.fasta > dates.txt

All taxa in the alignment are set to date 0 (present), except the three ancient
samples, which are set to -2677 (years before present, based on the radiocarbon
date of an associated faunal bone; see README.md). A calibration line for the
snowshoe hare-brown hare ancestral node (3.99 Mya) is written first.

Edit ANCIENT_SAMPLES, ANCIENT_DATE, CALIBRATION_TAXA and CALIBRATION_DATE below
if your sample names or calibration differ.
"""

import sys
from Bio import SeqIO

ANCIENT_SAMPLES = ["Lcor_Iron_125", "Lcor_Iron_127", "Lcor_Iron_129"]
ANCIENT_DATE = -2677  # years before present (775-541 BC calibrated radiocarbon
                       # midpoint, converted relative to contemporary sample
                       # collection date; see README.md and main text Methods)

# Taxa descending from the calibrated node (snowshoe hare + brown hare), and the
# age of their most recent common ancestor, in years before present.
CALIBRATION_TAXA = [
    "Lame1", "Lame2", "Lame3", "Lame4_pruned", "Lame5", "Lame6", "Lame7", "Lame8",
    "Leur2", "Leur3", "Leur4", "Leur5", "Leur6", "Leur7", "Leur8",
]
CALIBRATION_DATE = -3990000  # 3.99 Mya, Ferreira et al. 2021

def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: python3 make_dates_file.py alignment.fasta")

    names = [rec.id for rec in SeqIO.parse(sys.argv[1], "fasta")]
    if not names:
        sys.exit("No sequences found, check the input file")

    missing_calib = [t for t in CALIBRATION_TAXA if t not in names]
    if missing_calib:
        sys.exit(f"These calibration taxa are not in the alignment: {missing_calib}")
    missing_ancient = [t for t in ANCIENT_SAMPLES if t not in names]
    if missing_ancient:
        sys.exit(f"These ancient samples are not in the alignment: {missing_ancient}")

    print(",".join(CALIBRATION_TAXA) + f"\t{CALIBRATION_DATE}")
    for n in names:
        date = ANCIENT_DATE if n in ANCIENT_SAMPLES else 0
        print(f"{n}\t{date}")

if __name__ == "__main__":
    main()
