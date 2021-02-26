## PTC (premature termination codon) analysis of intron clusters
import os

import Intron_PTC_Annotation

# Annotate the PTC status of U intron
Intron_PTC_Annotation.get_intron_PTC()

# Calculate the proportion of PTC contained U introns for each intron clusters
Intron_PTC_Annotation.analyze_intron_cluster_ptc()

os.system('Rscript draw.ptc.R')
