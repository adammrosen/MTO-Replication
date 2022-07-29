
* This file creates all the Tables and Figures in Sections 5.1, 5.2 and Appendices B, C, D 
* plus CT estimates reported in Section 5.3 of "Estimating endogenous effects on ordinal outcomes"

* last updated: 	30/06/2022
* last executed: 	30/06/2022

* set options

clear all

set maxvar 10000
set linesize 255
set more off

set scheme s1color

* set directories for data files and output logs; change as needed

global pgmdir 	D:\data\12942365\ICPSR_34860\DS0001\
global datadir 	D:\data\12942365\ICPSR_34860\DS0001\data\
global pathres 	D:\data\12942365\ICPSR_34860\DS0001\stata_outputs\

* set date 

global date 30062022

* set dataset name

global data ${datadir}34860-0001-Data-REST

* open data set

use ${data}, clear

* recode female from male covariate

gen tmp_female = 1 - x_rad_ad_male

* define weight global

global wt f_wt_totsvy

* site covariates

global site_covs x_f_site_balt x_f_site_bos x_f_site_chi x_f_site_la

* set global "xcovs" with covariate variables (exogenous covariates)

global x_covariates 	x_rad_ad_ethrace_black_nh x_rad_ad_ethrace_hisp tmp_female x_rad_ad_le_35 x_rad_ad_36_40 x_rad_ad_41_45 x_rad_ad_46_50 ///
						x_f_ad_nevmarr x_f_ad_parentu18 x_f_ad_working x_f_ad_edinsch x_f_ad_edgradhs x_f_ad_edged ///
						x_f_hh_afdc x_f_hh_car x_f_hh_disabl x_f_hh_noteens x_f_hh_size2 x_f_hh_size3 x_f_hh_size4 x_f_hh_victim ///
						rad_svy_bl_totincm_2009d ///
						x_f_hood_unsafenit x_f_hood_verydissat x_f_hood_5y x_f_hood_nofamily x_f_hood_nofriend x_f_hood_chat x_f_hood_nbrkid ///
						x_f_hous_fndapt x_f_hous_sec8bef x_f_hous_mov3tm  x_f_hous_movdrgs x_f_hous_movschl cov_hous_movapt cov_hous_movjob 

* subjective well-being

gen happy_scale012 = happy_scale123_ad - 1

* generate the instrument, interaction between treatment group and site 

forvalues r = 1/3 {
	forvalues s = 1/5 {
		gen interaction_`r'_`s' = (ra_group==`r'&ra_site==`s')
	}
}

* set estimation samples

*regress happy_scale012 $site_covs f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z
*gen sample_nox = e(sample)

*regress happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z
*gen sample_x = e(sample)

gen sample_nox 	= 1
gen sample_x	= 1
	
* labels for randomization site and group 

label define site 1 "Baltimore" 2 "Boston" 3 "Chicago" 4 "LA" 5 "NY"
label values ra_site site

label define group 1 "Experimental" 2 "Section 8" 3 "Control"
label values ra_group group

** FIGURES 1-2: distribution of hood poverty and minority by randomization group 

do ${pgmdir}figures1-2.do

** APPENDIX B: descriptive statistics of baseline covariates by randomization group

do ${pgmdir}appendixB.do

** APPENDIX C: linear model estimates 

do ${pgmdir}appendixC.do

** TABLE 1: ordered probit estimates; requires myoprobit_lf.ado 

do ${pgmdir}table1.do

** table1_alt.do gives ordered logit and multinomial logit rather than ordered probit estimates, 
** these results are reported as a footnote in the paper. 
** This do-file also constructs figures showing probabilities and marginal effects using the alternative models  

do ${pgmdir}table1_alt.do

** APPENDIX D: CT estimates first stage; requires myoprobit2_lf.ado and myoprobit3_lf.ado 

do ${pgmdir}appendixD.do

** TABLE 2: CT estimates; requires myoprobit2_lf.ado and myoprobit3_lf.ado 

do ${pgmdir}table2.do

** table2_alt.do also gives the CT estimates but using a control function approach rather than maximum likelihood, 
** these results are reported as a footnote in the paper.  
** This do-file also constructs figures showing probabilities and marginal effects using the control function approach   

preserve 

do ${pgmdir}table2_alt.do

restore 
	
** FIGURES 3-6: constructs figures 3 to 6 in the paper; requires myoprobit_lf.ado, myoprobit2_lf.ado and myoprobit3_lf.ado 

preserve

do ${pgmdir}figures3-6.do

restore

** figures3-6_nox.do re-constructs figures 3 to 6 using specifications which do not control for baseline covariates, 
** these results are reported as a footnote in the paper.
** requires myoprobit_lf.ado, myoprobit2_lf.ado and myoprobit3_lf.ado 

preserve

do ${pgmdir}figures3-6_nox.do

restore

** InfoMatrixTests.do carries out the Information Matrix tests reported in the paper; requires myoprobit2_lf.ado and myoprobit3_lf.ado  

preserve

do ${pgmdir}InfoMatrixTests.do

restore

** CT in TABLES 3-5: estimates and reports the CT estimates in Tables 3-5; requires myoprobit2_lf.ado and myoprobit3_lf.ado 

preserve

do ${pgmdir}tables3-5-CT.do

restore

** end program 

disp "Program End on " c(current_date) " at " c(current_time)
