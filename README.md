# MTO-Replication
Replication Files for MTO Application of Chesher, Rosen, & Siddique (2023)

This repository contains replication code for the 2023 revision of "Estimating Endogenous Effects of Ordinal Outcomes" by Chesher, Rosen, and Siddique.

The U.S. Department of Housing and Urban Development (HUD) provided the MTO data; we used the version made available by the Inter-university Consortium for Political and Social Research (ICPSR) at the University of Michigan.

The data is not publicly available, but may be accessed with permission through ICPSR. See https://www.icpsr.umich.edu/web/ICPSR/studies/34860.
Researchers who obtain access to the data from ICPSR will need to set path variables accordingly to point to where the data is locally accessed. The local path used here was "D:/data/12942365/ICPSR_34860/DS0001".

The results presented for the numerical illustrations of Section 3.2 and Table 1 were obtained by executing the file "MTO_numerical_example_code.R" contained in the folder "Numerical Illustration".  Further details regarding the R code are contained in the README file in that folder.  These illustrations are based on hypothetical data generating processes and do not use MTO data.

For researchers who have access to the MTO data, code in the folder “RCPP” can be used to replicate set estimates and confidence sets reported in Tables 4–6.  Researchers who do not have access to the data may alternatively compute set estimates and confidence sets using simulated data.  See the comments at the beginning of “Profile_IV_MTO.R” for details.

Researchers with access to the data can replicate all other empirical results in the paper and supplementary materials with code in the folder “Stata”. All of these additional results can be replicated using "MTO_MAIN_final.do".
