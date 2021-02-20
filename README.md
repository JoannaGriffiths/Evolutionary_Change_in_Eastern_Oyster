# Scripts and data from Griffiths, J.S. et al. (in review) Experimental evolution reveals standing genetic variation for salinity tolerance in the eastern oyster.

## PCoA Analysis

*Script*: PCoA_oysters.R, treatments_all.txt, treatments_larvae.txt, PCoA_larval_distances.txt

*Input files*: all_pop_exact_cov20_200_new_rc_edit2

*Description*: Script contains code for running a Principal Coordinate Analysis for all samples (parents and larvae) and larvae samples only. The script also contains adonis functions for testing significance of groups. The treatments file contains treatment info for all samples. The distance file contains Manhattan distances for larval samples pre- and post-low salinity selection from the PCoA plot.




## Detecting increases in allele frequencies

*Script*: Manhattan_plot.R 

*Input files*: *CROSS*-newpvalue_Manhattan (e.g. LA1-newpvalue_Manhattan)

*Description*: Script contains code for plotting Manhattan plots for crosses that had genes under selection and crosses that did not.




## Annotation of Genes Under Selection

*Script*: Annotation.R

*Input files*: 

*Description*: Script contains code for determining if genes under seleciton were enriched in a particular functional or source category. It also determines whether SNPs were more likely to be upstream of the gene than within the gene body.

