# Import modules
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import diempy as diem
import argparse
import subprocess
import os

# Parse arguments
parser = argparse.ArgumentParser(description="Run DIEM pipeline from VCF")
parser.add_argument("vcf", help="Input VCF file (.vcf.gz)")
parser.add_argument("-c", "--cores", type=int, default=10,
                    help="Number of CPU cores to use (default: 10)")
args = parser.parse_args()

vcfFile = args.vcf
ncores = args.cores
prefix = vcfFile

# Run vcf2diem if diem input doesn't exist
bedInFile = prefix + ".diem_input.bed"
metaFile = prefix + ".diem_meta.bed"
if not os.path.exists(bedInFile):
    subprocess.run(["vcf2diem", vcfFile], check=True)

# Read DIEM BED files
dRaw = diem.read_diem_bed(bedInFile, metaFile)
print("Successfully loaded:", prefix)

# Polarize if polarized file doesn't exist
polarizedDiemTypeOutFile = prefix + '.polarized.diemtype'
if os.path.exists(polarizedDiemTypeOutFile):
    dPol = diem.load_DiemType(polarizedDiemTypeOutFile)
    print(f"Loaded polarized DIEM object from {polarizedDiemTypeOutFile}")
else:
    dPol = dRaw.polarize(ncores=ncores)
    print(f"Polarization complete using {ncores} cores.")
    # Save polarized objects
    diem.write_polarized_bed(bedInFile, prefix + '.diem_output.bed', dPol)
    diem.save_DiemType(dPol, polarizedDiemTypeOutFile)

# Sort individuals
dPol.sort()

# Plot HI for raw (no threshold)
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(dPol.HIs, marker='.')
ax.set_ylim(0, 1)
ax.set_ylabel('HI')
ax.set_xlabel('Individual')
ax.set_xticks(range(len(dPol.indNames)))
ax.set_xticklabels(dPol.indNames, rotation=90)
plt.savefig(prefix + "_HI_plot.pdf", bbox_inches="tight")
plt.close(fig)

# Save raw HI
dfHI = pd.DataFrame(list(zip(dPol.indNames, dPol.HIs)),
                    columns=["Individual", "HI"])
dfHI.to_csv(f"{prefix}_HI_per_individual_raw.tsv", sep="\t", index=False)
print(f"Saved HI table to: {prefix}_HI_per_individual_raw.tsv")

# Thresholds to apply
thresholds = [-40, -20, -15, -10, None]

for thresh in thresholds:
    if thresh is not None:
        dThresh = dPol.apply_threshold(thresh)
        hi_vals = list(zip(dPol.indNames, dThresh.HIs))
        tsv_outfile = f"{prefix}_HI_after_threshold_{thresh}.tsv"
    else:
        hi_vals = list(zip(dPol.indNames, dPol.HIs))
        tsv_outfile = f"{prefix}_HI_no_threshold.tsv"
    
    dfHI_thresh = pd.DataFrame(hi_vals, columns=["Individual", "HI"])
    dfHI_thresh.to_csv(tsv_outfile, sep="\t", index=False)
    print(f"Saved HI table for threshold {thresh} to: {tsv_outfile}")
