# MTO-Replication
 Replication Files for MTO Application of Chesher, Rosen, & Siddique (2022)

This repository contains replication code for the 2022 working paper "Estimating Endogenous Effects of Ordinal Outcomes" by Chesher, Rosen, and Siddique.

The U.S. Department of Housing and Urban Development (HUD) provided the MTO data; we used the version made available by the Inter-university Consortium for Political and Social Research (ICPSR) at the University of Michigan.

The data is not publicly available, but may be accessed with permission through ICPSR. See https://www.icpsr.umich.edu/web/ICPSR/studies/34860.
Researchers who obtain access to the data from ICPSR will need to set path variables accordingly to point to where the data is locally accessed. The local path used here was "D:/data/12942365/ICPSR_34860/DS0001/R".

For researchers who have access to the data, code in the folder “RCPP” can be used to replicate set estimates and confidence sets reported in Tables 3–5.  Researchers who do not have access to the data may alternatively compute set estimates and confidence sets using simulated data.  See the comments at the beginning of “Profile_IV_MTO.R” for details.

Researchers with access to the data can replicate all other results in the paper and supplementary materials with code in the folder “Stata”. All results can be replicated using "MTO_MAIN_final.do".
