Bayesian tip-dating (BEAST2, relaxed lognormal clock, node + tip calibrations combined)
gave unstable rate estimates when the ancient samples were added, because it had to
reconcile calibrations of very different scales (a node in the millions of years, tips a
few thousand years old). This is a known issue, usually called time-dependent rate bias
(Ho et al. 2005; Ho et al. 2011). Switched to least-squares dating (LSD2, via IQ-TREE3) for
the main analysis, see 02.IQTREE_Dating.sh. The original BEAST2 analysis, restricted to
the contemporary samples only, is kept as a confirmatory analysis (Supplementary Text S1,
mito_hares.xml in this folder).
