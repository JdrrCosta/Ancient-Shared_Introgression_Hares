#Import modules
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
np.set_printoptions(legacy = '1.25')

import diempy as diem
import time

import sys
import os

import argparse
import subprocess

#parse vcf
parser = argparse.ArgumentParser(description="Run DIEM pipeline from VCF")
parser.add_argument("vcf", help="Input VCF file (.vcf.gz)")
parser.add_argument("-c", "--cores", type=int, default=10,
                    help="Number of CPU cores to use (default: 10)")
                    
args = parser.parse_args()

vcfFile = args.vcf
ncores = args.cores

#run vcf2diem
subprocess.run(["vcf2diem", vcfFile], check=True)

#derive prefixes
prefix = vcfFile

bedInFile = prefix + ".diem_input.bed"
metaFile = prefix + ".diem_meta.bed"

#mask individuals
individualsMaskedFile = './individuals_masked.txt'
#read into diempy
dRawWithMasks = diem.read_diem_bed(bedInFile, metaFile)
dRawWithMasks.add_individual_exclusions(individualsMaskedFile)

print("Successfully loaded:", prefix)

#polarization
dPol = dRawWithMasks.polarize(ncores=ncores)
print(f"Polarization complete using {ncores} cores.")

#save to bed file
bedOutFile = prefix + '.diem_output.bed'
diem.write_polarized_bed(bedInFile,bedOutFile,dPol)
#save this as a diemtype object
polarizedDiemTypeOutFile = prefix + '.polarized.diemtype'
diem.save_DiemType(dPol,polarizedDiemTypeOutFile)

#sort indv by hibdrid index
dPol.sort()
#plot hybrid indexes
fig, ax = plt.subplots(figsize=(6,4))
ax.plot(dPol.HIs, marker='.')
ax.set_ylim(0, 1)
ax.set_ylabel('HI')
ax.set_xlabel('individual')

# Use individual names instead of numeric index
ax.set_xticks(range(len(dPol.indNames)))
ax.set_xticklabels(dPol.indNames, rotation=90)  # rotate for readability

plt.savefig(prefix + "_HI_plot.pdf", bbox_inches="tight")
plt.close(fig)

#list of HI before filtering
#Create a table of HI per individual
hi_per_ind = list(zip(dPol.indNames, dPol.HIs))
# Convert to a pandas DataFrame
dfHI = pd.DataFrame(hi_per_ind, columns=["Individual", "HI"])
#Save as TSV for downstream analysis
out_file = f"{prefix}_HI_per_individual_raw.tsv"
dfHI.to_csv(out_file, sep="\t", index=False)

print(f"Saved HI table to: {out_file}")

#plot diem painting
#here we are plotting only one chromossome
diem.plot_painting(dPol.DMBC[0],names=dPol.indNames)

#plot diagnostic index DI
DIValsWithRepeats = np.hstack(dPol.DIByChr)
DIValsUnique = np.unique(DIValsWithRepeats)

fig, ax = plt.subplots(figsize=(5,4))

counts, bins, patches = ax.hist(
    DIValsWithRepeats,
    bins=500,
    color='darkmagenta'
)

ax.hist(
    DIValsUnique,
    bins=bins,
    color='lightseagreen'
)

ax.set_yscale('log')

fig.savefig(prefix + "_DI_histogram.pdf", bbox_inches="tight")
plt.close(fig)

#charactrize markers with a threshold above -50 (could be computer intensive)
dfMarkers = diem.characterize_markers(dPol.apply_threshold(-50))
dfMarkers.to_csv(prefix + "_dfMarkers.csv", index=False)

#set a DI distribution treshold
dThresh = dPol.apply_threshold(-20)

#plot before and after thresholding
# ---- BEFORE threshold ----
diem.plot_painting_with_positions(
    dPol.DMBC[0],
    dPol.posByChr[0],
    names=dPol.indNames
)

plt.savefig(prefix + "_painting_before_threshold.pdf", bbox_inches="tight")
plt.close()


# ---- AFTER threshold ----
diem.plot_painting_with_positions(
    dThresh.DMBC[0],
    dThresh.posByChr[0],
    names=dPol.indNames
)

plt.savefig(prefix + "_painting_after_threshold.pdf", bbox_inches="tight")
plt.close()

# Hybrid index before and after thresholding
fig, ax = plt.subplots(figsize=(6,4))

ax.plot(dPol.HIs, marker='.', color='b', label='Before threshold')
ax.plot(dThresh.HIs, marker='.', color='r', label='After threshold')

ax.set_ylim(0, 1)
ax.set_ylabel('HI')
ax.set_xlabel('individual')
ax.legend()

# Use individual names instead of numeric index
ax.set_xticks(range(len(dPol.indNames)))
ax.set_xticklabels(dPol.indNames, rotation=90)

plt.tight_layout()  # avoid label clipping
fig.savefig(prefix + "_HI_comparison.pdf", bbox_inches="tight")
plt.close(fig)

#sve HI after thresholding
hi_after_thresh = list(zip(dPol.indNames, dThresh.HIs))
dfHI_thresh = pd.DataFrame(hi_after_thresh, columns=["Individual", "HI_after_threshold"])

tsv_outfile = f"{prefix}_HI_after_threshold.tsv"
dfHI_thresh.to_csv(tsv_outfile, sep="\t", index=False)

print(f"Saved HI after thresholding table to: {tsv_outfile}")

#smoothing the plot
scalesToTry = [1e-8,5e-8,1e-7,5e-7,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4]
sitesDiffByScale = []
for scale in scalesToTry:
    dSmoothedTest = dThresh.smooth(scale)
    kdiffs = diem.count_site_differences(dSmoothedTest.DMBC,dThresh.DMBC)
    sitesDiffByScale.append(kdiffs)
    
#number of sites changed
# Create the figure
fig, ax = plt.subplots(figsize=(6,4))

# Plot the data
ax.plot(scalesToTry, sitesDiffByScale, marker='.', linestyle='')  # markers only

# Set axis labels
ax.set_ylabel('# of sites changed')
ax.set_xlabel('Laplace scale')

# Adjust layout to prevent clipping
plt.tight_layout()

# Save as PDF
fig.savefig(prefix + "_sites_diff_by_scale.pdf", bbox_inches="tight")

# Close the figure to free memory
plt.close(fig)


#pick smoothing threshold
dSmoothedUnsorted = dThresh.smooth(1e-5)

#compare plots before and after smoothing
# ---- Plot after thresholding ----
diem.plot_painting(dThresh.DMBC[0], dThresh.indNames)

plt.tight_layout()  # ensure layout is nice
plt.savefig(prefix + "_painting_after_threshold.pdf", bbox_inches="tight")
plt.close()  # free memory

# ---- Plot smoothed unsorted ----
diem.plot_painting(dSmoothedUnsorted.DMBC[0], dSmoothedUnsorted.indNames)

plt.tight_layout()
plt.savefig(prefix + "_painting_smoothed_unsorted.pdf", bbox_inches="tight")
plt.close()

#HI before and after smoothing
# Create figure
fig, ax = plt.subplots(figsize=(6,4))

# Plot HI values
ax.plot(dThresh.HIs, marker='.', color='b', label='After threshold')
ax.plot(dSmoothedUnsorted.computeHIs(), marker='.', color='r', label='Smoothed unsorted')

# Set labels
ax.set_ylabel('HI')
ax.set_xlabel('individual idx')
ax.legend()

# Adjust layout to avoid clipping
plt.tight_layout()

# Save figure to PDF
fig.savefig(prefix + "_HI_comparison_thresh_vs_smoothed.pdf", bbox_inches="tight")

# Close figure to free memory
plt.close(fig)

hi_smoothed = dSmoothedUnsorted.computeHIs()
hi_table = list(zip(dPol.indNames, dThresh.HIs, hi_smoothed))
dfHI = pd.DataFrame(hi_table, columns=["Individual", "HI_after_threshold", "HI_smoothed"])

tsv_outfile = f"{prefix}_HI_thresh_vs_smoothed.tsv"
dfHI.to_csv(tsv_outfile, sep="\t", index=False)

print(f"Saved HI table for thresholded vs smoothed values to: {tsv_outfile}")

#smooothed sort
# Make a sorted copy
dSmoothed = dSmoothedUnsorted.copy()
dSmoothed.sort()

# Plot the painting
diem.plot_painting(dSmoothed.DMBC[0], dSmoothed.indNames)

# Adjust layout and save as PDF
plt.tight_layout()
plt.savefig(prefix + "_painting_smoothed_sorted.pdf", bbox_inches="tight")

# Close figure to free memory
plt.close()


# Create the DeFinetti plot
fig = diem.GenomicDeFinettiPlot(dPol)

# Save the current active figure
plt.savefig(prefix + "_DeFinetti_plot.pdf", bbox_inches="tight")

# Close it to free memory
plt.close()

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

#DeFinetti plots (chr1–23)
chroms_to_plot = list(range(23))  # 0–20 = chr1–23

diem.GenomicMultiDeFinettiPlot(dPol, chroms_to_plot)
plt.gcf().savefig(prefix + "_GenomicMultiDeFinetti_chr1-23.pdf", bbox_inches="tight")
plt.close()


# -----------------------------
# Genomic contributions (chr1–23)
# -----------------------------
diem.GenomicContributionsPlot(dPol, chroms_to_plot)
plt.gcf().savefig(prefix + "_GenomicContributions_chr1-23.pdf", bbox_inches="tight")
plt.close()


# -----------------------------
# Individual genomic contributions
# -----------------------------
diem.IndGenomicContributionsPlot(dPol)
plt.gcf().savefig(prefix + "_IndGenomicContributions.pdf", bbox_inches="tight")
plt.close()


# -----------------------------
# PlotPrep object (ONLY define once)
# -----------------------------
MyFirstPlotPrep = diem.diemPlotPrepFromBedMeta(
    plot_theme="tutorial example",
    bed_file_path=bedOutFile,
    meta_file_path=metaFile,
    di_threshold=-57,
    genome_pixels=360*10,
    ticks=10000,
    smooth=10e-5
)


# -----------------------------
# Iris plot
# -----------------------------
diem.diemIrisFromPlotPrep(MyFirstPlotPrep)
plt.gcf().savefig(prefix + "_Iris_plot.pdf", bbox_inches="tight")
plt.close()


# -----------------------------
# Long plots chr1–23
# -----------------------------
for i in range(23):
    diem.diemLongFromPlotPrep(MyFirstPlotPrep, [i])
    plt.gcf().savefig(f"{prefix}_Long_chr{i+1}.pdf", bbox_inches="tight")
    plt.close()
