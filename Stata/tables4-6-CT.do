
** set scalars and globals to use later

sum f_c9010t_perpov_dw_z if (ra_site==5 & sample_nox==1), detail
scalar f_c9010t_perpov_dw_z_med_nox 	= r(p50)
scalar f_c9010t_perpov_dw_z_sd_nox 		= r(sd)

sum f_c9010t_pminorty_dw_z if (ra_site==5 & sample_nox==1), detail
scalar f_c9010t_pminorty_dw_z_med_nox 	= r(p50)
scalar f_c9010t_pminorty_dw_z_sd_nox 	= r(sd)

sum f_c9010t_perpov_dw_z if (ra_site==5 & sample_x==1), detail
scalar f_c9010t_perpov_dw_z_med_x 		= r(p50)
scalar f_c9010t_perpov_dw_z_sd_x 		= r(sd)

sum f_c9010t_pminorty_dw_z if (ra_site==5 & sample_x==1), detail
scalar f_c9010t_pminorty_dw_z_med_x 	= r(p50)
scalar f_c9010t_pminorty_dw_z_sd_x 		= r(sd)

foreach x in $x_covariates {		
	summarize `x' if (ra_site==5 & sample_x==1), detail
	scalar `x'_med_x = r(p50)
}
	
sum x_f_release1 if (ra_site==5 & sample_nox==1), detail
scalar x_f_release1_med_nox 			= r(p50)

sum x_f_release1 if (ra_site==5 & sample_x==1), detail
scalar x_f_release1_med_x 				= r(p50)

scalar x_rad_ad_41_45_med_x 			= 1
scalar x_f_hh_size3_med_x 				= 1

scalar x_rad_ad_ethrace_black_nh_med_x 	= 0 
scalar x_rad_ad_ethrace_hisp_med_x 		= 1 

global xgamma_short "(0) "
global xgamma_morex "x_f_release1_med_x*[xb1]x_f_release1+tmp_female_med_x*[xb1]tmp_female+x_rad_ad_le_35_med_x*[xb1]x_rad_ad_le_35+x_rad_ad_36_40_med_x*[xb1]x_rad_ad_36_40+x_rad_ad_41_45_med_x*[xb1]x_rad_ad_41_45+x_rad_ad_46_50_med_x*[xb1]x_rad_ad_46_50+x_rad_ad_ethrace_black_nh_med_x*[xb1]x_rad_ad_ethrace_black_nh+x_rad_ad_ethrace_hisp_med_x*[xb1]x_rad_ad_ethrace_hisp+x_f_ad_nevmarr_med_x*[xb1]x_f_ad_nevmarr+x_f_ad_parentu18_med_x*[xb1]x_f_ad_parentu18+x_f_ad_working_med_x*[xb1]x_f_ad_working+x_f_ad_edinsch_med_x*[xb1]x_f_ad_edinsch+x_f_ad_edgradhs_med_x*[xb1]x_f_ad_edgradhs+x_f_ad_edged_med_x*[xb1]x_f_ad_edged+x_f_hh_afdc_med_x*[xb1]x_f_hh_afdc+rad_svy_bl_totincm_2009d_med_x*[xb1]rad_svy_bl_totincm_2009d+x_f_hh_car_med_x*[xb1]x_f_hh_car+x_f_hh_disabl_med_x*[xb1]x_f_hh_disabl+x_f_hh_noteens_med_x*[xb1]x_f_hh_noteens+x_f_hh_size2_med_x*[xb1]x_f_hh_size2+x_f_hh_size3_med_x*[xb1]x_f_hh_size3+x_f_hh_size4_med_x*[xb1]x_f_hh_size4+x_f_hh_victim_med_x*[xb1]x_f_hh_victim+x_f_hood_unsafenit_med_x*[xb1]x_f_hood_unsafenit+x_f_hood_verydissat_med_x*[xb1]x_f_hood_verydissat+x_f_hood_5y_med_x*[xb1]x_f_hood_5y+x_f_hous_mov3tm_med_x*[xb1]x_f_hous_mov3tm+x_f_hood_nofamily_med_x*[xb1]x_f_hood_nofamily+x_f_hood_nofriend_med_x*[xb1]x_f_hood_nofriend+x_f_hood_chat_med_x*[xb1]x_f_hood_chat+x_f_hood_nbrkid_med_x*[xb1]x_f_hood_nbrkid+x_f_hous_fndapt_med_x*[xb1]x_f_hous_fndapt+x_f_hous_sec8bef_med_x*[xb1]x_f_hous_sec8bef+x_f_hous_movdrgs_med_x*[xb1]x_f_hous_movdrgs+x_f_hous_movschl_med_x*[xb1]x_f_hous_movschl+cov_hous_movapt_med_x*[xb1]cov_hous_movapt+cov_hous_movjob_med_x*[xb1]cov_hous_movjob"

** set value of hood poverty to median in order to estimate the CT marginal effects reported in Tables 3 and 4
** as well as counterfactual probabilities reported in columns 1 and 2 of Table 5

scalar hood_pov_touse_nox = f_c9010t_perpov_dw_z_med_nox
scalar hood_pov_touse_x = f_c9010t_perpov_dw_z_med_x

** TABLE 4

file open csvlog using "${pathres}table4_${date}.csv", write replace

** CT model (W=hood poverty) without covariates 

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs f_c9010t_perpov_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_perpov_dw_z = interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml init xb1:x_f_site_bal=0.1825 x_f_site_bos=0.0528 x_f_site_chi=0.1811 x_f_site_la=0.0674 f_c9010t_perpov_dw_z=-0.1610
ml maximize, difficult

*** marginal effects for Y = 0,1,2 (reported in column 1 of Table 4)

file write csvlog "y, ME(Y=y|pov), lower CI, upper CI" _n

nlcom -[xb1]f_c9010t_perpov_dw_z * normalden([cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}))
matrix b = r(b)
matrix V = r(V)
scalar me_y0_povmed = b[1,1]
scalar me_y0_povmed_v = sqrt(V[1,1])
scalar me_y0_povmed_lci = me_y0_povmed-invnormal(0.975)*me_y0_povmed_v
scalar me_y0_povmed_uci = me_y0_povmed+invnormal(0.975)*me_y0_povmed_v

disp "marginal effect that y=0 is " me_y0_povmed " and the 95% confidence interval is from "me_y0_povmed_lci " to " me_y0_povmed_uci
file write csvlog " 0 ," %20.5f (me_y0_povmed) "," %20.5f (me_y0_povmed_lci) "," %20.5f (me_y0_povmed_uci) "," _n

nlcom [xb1]f_c9010t_perpov_dw_z * (normalden([cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}))- /*
				*/normalden(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short})))
matrix b = r(b)
matrix V = r(V)
scalar me_y1_povmed = b[1,1]
scalar me_y1_povmed_v = sqrt(V[1,1])
scalar me_y1_povmed_lci = me_y1_povmed-invnormal(0.975)*me_y1_povmed_v
scalar me_y1_povmed_uci = me_y1_povmed+invnormal(0.975)*me_y1_povmed_v

disp "marginal effect that y=1 is " me_y1_povmed " and the 95% confidence interval is from "me_y1_povmed_lci " to " me_y1_povmed_uci
file write csvlog " 1 ," %20.5f (me_y1_povmed) "," %20.5f (me_y1_povmed_lci) "," %20.5f (me_y1_povmed_uci) "," _n

nlcom [xb1]f_c9010t_perpov_dw_z * normalden(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short})) 
matrix b = r(b)
matrix V = r(V)
scalar me_y2_povmed = b[1,1]
scalar me_y2_povmed_v = sqrt(V[1,1])
scalar me_y2_povmed_lci = me_y2_povmed-invnormal(0.975)*me_y2_povmed_v
scalar me_y2_povmed_uci = me_y2_povmed+invnormal(0.975)*me_y2_povmed_v

disp "marginal effect that y=2 is " me_y2_povmed " and the 95% confidence interval is from "me_y2_povmed_lci " to " me_y2_povmed_uci
file write csvlog " 2 ," %20.5f (me_y2_povmed) "," %20.5f (me_y2_povmed_lci) "," %20.5f (me_y2_povmed_uci) "," _n

** CT model (W= hood poverty) with covariates 

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_perpov_dw_z = x_f_release1 $x_covariates interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml init xb1:x_f_site_bal=0.1825 x_f_site_bos=0.0528 x_f_site_chi=0.1811 x_f_site_la=0.0674 f_c9010t_perpov_dw_z=-0.1610
ml maximize, difficult

*** marginal effects for Y = 0,1,2 (reported in column 2 of Table 4)

file write csvlog "y, ME(Y=y|povX), lower CI, upper CI" _n

nlcom -[xb1]f_c9010t_perpov_dw_z * normalden([cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}+ ${xgamma_morex}) )
matrix b = r(b)
matrix V = r(V)
scalar me_x_y0_povmed = b[1,1]
scalar me_x_y0_povmed_v = sqrt(V[1,1])
scalar me_x_y0_povmed_lci = me_x_y0_povmed-invnormal(0.975)*me_x_y0_povmed_v
scalar me_x_y0_povmed_uci = me_x_y0_povmed+invnormal(0.975)*me_x_y0_povmed_v

disp "marginal effect that y=0 is " me_x_y0_povmed " and the 95% confidence interval is from "me_x_y0_povmed_lci " to " me_x_y0_povmed_uci
file write csvlog " 0 ," %20.5f (me_x_y0_povmed) "," %20.5f (me_x_y0_povmed_lci) "," %20.5f (me_x_y0_povmed_uci) "," _n

nlcom [xb1]f_c9010t_perpov_dw_z * (normalden([cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}+ ${xgamma_morex}))- /*
				*/normalden(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}+ ${xgamma_morex}) ))
matrix b = r(b)
matrix V = r(V)
scalar me_x_y1_povmed = b[1,1]
scalar me_x_y1_povmed_v = sqrt(V[1,1])
scalar me_x_y1_povmed_lci = me_x_y1_povmed-invnormal(0.975)*me_x_y1_povmed_v
scalar me_x_y1_povmed_uci = me_x_y1_povmed+invnormal(0.975)*me_x_y1_povmed_v

disp "marginal effect that y=1 is " me_y1_povmed " and the 95% confidence interval is from "me_y1_povmed_lci " to " me_y1_povmed_uci
file write csvlog " 1 ," %20.5f (me_x_y1_povmed) "," %20.5f (me_x_y1_povmed_lci) "," %20.5f (me_x_y1_povmed_uci) "," _n

nlcom [xb1]f_c9010t_perpov_dw_z * normalden(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}+ ${xgamma_morex}) )  
matrix b = r(b)
matrix V = r(V)
scalar me_x_y2_povmed = b[1,1]
scalar me_x_y2_povmed_v = sqrt(V[1,1])
scalar me_x_y2_povmed_lci = me_x_y2_povmed-invnormal(0.975)*me_x_y2_povmed_v
scalar me_x_y2_povmed_uci = me_x_y2_povmed+invnormal(0.975)*me_x_y2_povmed_v

disp "marginal effect that y=2 is " me_y2_povmed " and the 95% confidence interval is from "me_y2_povmed_lci " to " me_y2_povmed_uci
file write csvlog " 2 ," %20.5f (me_x_y2_povmed) "," %20.5f (me_x_y2_povmed_lci) "," %20.5f (me_x_y2_povmed_uci) "," _n

file close csvlog

** TABLE 5

file open csvlog using "${pathres}table5_${date}.csv", write replace

** CT model (W= hood poverty, hood minority) without covariates 

ml model lf myoprobit3_lf (xb1:happy_scale012 = $site_covs f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z, nocons) (cut1:) (lndiff:) /*
						*/(xb2:f_c9010t_perpov_dw_z = interaction*) (lnvar_v1:) (rho_uv1:) /*
						*/(xb3:f_c9010t_pminorty_dw_z = interaction*) (lnvar_v2:) (rho_uv2:) (cov_v1v2:) [pw=$wt], /*
						*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
						*/diparm(lnvar_v1, exp label(var_v1)) /*
						*/diparm(lnvar_v2, exp label(var_v2)) 
ml init xb1:x_f_site_bal=0.3039 x_f_site_bos=0.3386 x_f_site_chi=0.1527 x_f_site_la=0.0639 f_c9010t_perpov_dw_z=-0.2998 f_c9010t_pminorty_dw_z=0.3151 cut1:_cons=-0.8165
ml maximize, difficult

*** marginal effects for Y = 0,1,2 (reported in column 1 of Table 5)

file write csvlog "y, ME(Y=y|povmin), lower CI, upper CI" _n

nlcom -[xb1]f_c9010t_perpov_dw_z * normalden([cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_nox*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}))
matrix b = r(b)
matrix V = r(V)
scalar me_y0_povminmed = b[1,1]
scalar me_y0_povminmed_v = sqrt(V[1,1])
scalar me_y0_povminmed_lci = me_y0_povminmed-invnormal(0.975)*me_y0_povminmed_v
scalar me_y0_povminmed_uci = me_y0_povminmed+invnormal(0.975)*me_y0_povminmed_v

disp "marginal effect that y=0 is " me_y0_povminmed " and the 95% confidence interval is from "me_y0_povminmed_lci " to " me_y0_povminmed_uci
file write csvlog " 0 ," %20.5f (me_y0_povminmed) "," %20.5f (me_y0_povminmed_lci) "," %20.5f (me_y0_povminmed_uci) "," _n

nlcom [xb1]f_c9010t_perpov_dw_z * (normalden([cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_nox*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}))- /*
				*/normalden(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_nox*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short})))
matrix b = r(b)
matrix V = r(V)
scalar me_y1_povminmed = b[1,1]
scalar me_y1_povminmed_v = sqrt(V[1,1])
scalar me_y1_povminmed_lci = me_y1_povminmed-invnormal(0.975)*me_y1_povminmed_v
scalar me_y1_povminmed_uci = me_y1_povminmed+invnormal(0.975)*me_y1_povminmed_v

disp "marginal effect that y=1 is " me_y1_povminmed " and the 95% confidence interval is from "me_y1_povminmed_lci " to " me_y1_povminmed_uci
file write csvlog " 1 ," %20.5f (me_y1_povminmed) "," %20.5f (me_y1_povminmed_lci) "," %20.5f (me_y1_povminmed_uci) "," _n

nlcom [xb1]f_c9010t_perpov_dw_z * normalden(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_nox*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}))
matrix b = r(b)
matrix V = r(V)
scalar me_y2_povminmed = b[1,1]
scalar me_y2_povminmed_v = sqrt(V[1,1])
scalar me_y2_povminmed_lci = me_y2_povminmed-invnormal(0.975)*me_y2_povminmed_v
scalar me_y2_povminmed_uci = me_y2_povminmed+invnormal(0.975)*me_y2_povminmed_v

disp "marginal effect that y=2 is " me_y2_povminmed " and the 95% confidence interval is from "me_y2_povminmed_lci " to " me_y2_povminmed_uci
file write csvlog " 2 ," %20.5f (me_y2_povminmed) "," %20.5f (me_y2_povminmed_lci) "," %20.5f (me_y2_povminmed_uci) "," _n

** CT model (W= hood poverty, hood minority) with covariates 
 
ml model lf myoprobit3_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z, nocons) (cut1:) (lndiff:) /*
						*/(xb2:f_c9010t_perpov_dw_z = x_f_release1 $x_covariates interaction*) (lnvar_v1:) (rho_uv1:) /*
						*/(xb3:f_c9010t_pminorty_dw_z = x_f_release1 $x_covariates interaction*) (lnvar_v2:) (rho_uv2:) (cov_v1v2:) [pw=$wt], /*
						*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
						*/diparm(lnvar_v1, exp label(var_v1)) /*
						*/diparm(lnvar_v2, exp label(var_v2)) 
ml init xb1:x_f_site_bal=0.3039 x_f_site_bos=0.3386 x_f_site_chi=0.1527 x_f_site_la=0.0639 f_c9010t_perpov_dw_z=-0.2998 f_c9010t_pminorty_dw_z=0.3151 cut1:_cons=-0.8165
ml maximize, difficult
est store povmin_t_x

*** marginal effects for Y = 0,1,2 (reported in column 2 of Table 5)

file write csvlog "y, ME(Y=y|povminX), lower CI, upper CI" _n

nlcom -[xb1]f_c9010t_perpov_dw_z * normalden([cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_x*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))
matrix b = r(b)
matrix V = r(V)
scalar me_x_y0_povminmed = b[1,1]
scalar me_x_y0_povminmed_v = sqrt(V[1,1])
scalar me_x_y0_povminmed_lci = me_x_y0_povminmed-invnormal(0.975)*me_x_y0_povminmed_v
scalar me_x_y0_povminmed_uci = me_x_y0_povminmed+invnormal(0.975)*me_x_y0_povminmed_v

disp "marginal effect that y=0 is " me_x_y0_povminmed " and the 95% confidence interval is from "me_x_y0_povminmed_lci " to " me_x_y0_povminmed_uci
file write csvlog " 0 ," %20.5f (me_x_y0_povminmed) "," %20.5f (me_x_y0_povminmed_lci) "," %20.5f (me_x_y0_povminmed_uci) "," _n

nlcom [xb1]f_c9010t_perpov_dw_z * normalden([cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_x*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))- /*
				*/[xb1]f_c9010t_perpov_dw_z * normalden(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_x*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))
matrix b = r(b)
matrix V = r(V)
scalar me_x_y1_povminmed = b[1,1]
scalar me_x_y1_povminmed_v = sqrt(V[1,1])
scalar me_x_y1_povminmed_lci = me_x_y1_povminmed-invnormal(0.975)*me_x_y1_povminmed_v
scalar me_x_y1_povminmed_uci = me_x_y1_povminmed+invnormal(0.975)*me_x_y1_povminmed_v

disp "marginal effect that y=1 is " me_x_y1_povminmed " and the 95% confidence interval is from "me_x_y1_povminmed_lci " to " me_x_y1_povminmed_uci
file write csvlog " 1 ," %20.5f (me_x_y1_povminmed) "," %20.5f (me_x_y1_povminmed_lci) "," %20.5f (me_x_y1_povminmed_uci) "," _n

nlcom [xb1]f_c9010t_perpov_dw_z * normalden(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_x*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))
matrix b = r(b)
matrix V = r(V)
scalar me_x_y2_povminmed = b[1,1]
scalar me_x_y2_povminmed_v = sqrt(V[1,1])
scalar me_x_y2_povminmed_lci = me_x_y2_povminmed-invnormal(0.975)*me_x_y2_povminmed_v
scalar me_x_y2_povminmed_uci = me_x_y2_povminmed+invnormal(0.975)*me_x_y2_povminmed_v

disp "marginal effect that y=2 is " me_x_y2_povminmed " and the 95% confidence interval is from "me_x_y2_povminmed_lci " to " me_x_y2_povminmed_uci
file write csvlog " 2 ," %20.5f (me_x_y2_povminmed) "," %20.5f (me_x_y2_povminmed_lci) "," %20.5f (me_x_y2_povminmed_uci) "," _n

file close csvlog

** TABLE 6a: columns 1 and 2 

file open csvlog using "${pathres}table6a_${date}.csv", write replace

** CT model (W=hood poverty) without covariates 

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs f_c9010t_perpov_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_perpov_dw_z = interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml init xb1:x_f_site_bal=0.1825 x_f_site_bos=0.0528 x_f_site_chi=0.1811 x_f_site_la=0.0674 f_c9010t_perpov_dw_z=-0.1610
ml maximize, difficult

*** counterfactual probabilities that Y = 0,1,2 (reported in column 1 of Table 6)

file write csvlog "y, P(Y=y|pov), lower CI, upper CI" _n

nlcom normal([cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}))
matrix b = r(b)
matrix V = r(V)
scalar pr_y0_povmed = b[1,1]
scalar pr_y0_povmed_v = sqrt(V[1,1])
scalar pr_y0_povmed_lci = pr_y0_povmed-invnormal(0.975)*pr_y0_povmed_v
scalar pr_y0_povmed_uci = pr_y0_povmed+invnormal(0.975)*pr_y0_povmed_v

disp "probability that y=0 is " pr_y0_povmed " and the 95% confidence interval is from "pr_y0_povmed_lci " to " pr_y0_povmed_uci
file write csvlog " 0 ," %20.5f (pr_y0_povmed) "," %20.5f (pr_y0_povmed_lci) "," %20.5f (pr_y0_povmed_uci) "," _n

nlcom normal(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short})) - /*
							*/normal([cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}))
matrix b = r(b)
matrix V = r(V)
scalar pr_y1_povmed = b[1,1]
scalar pr_y1_povmed_v = sqrt(V[1,1])
scalar pr_y1_povmed_lci = pr_y1_povmed-invnormal(0.975)*pr_y1_povmed_v
scalar pr_y1_povmed_uci = pr_y1_povmed+invnormal(0.975)*pr_y1_povmed_v

disp "probability that y=1 is " pr_y1_povmed " and the 95% confidence interval is from "pr_y1_povmed_lci " to " pr_y1_povmed_uci
file write csvlog " 1 ," %20.5f (pr_y1_povmed) "," %20.5f (pr_y1_povmed_lci) "," %20.5f (pr_y1_povmed_uci) "," _n

nlcom 1-normal(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short})) 
matrix b = r(b)
matrix V = r(V)
scalar pr_y2_povmed = b[1,1]
scalar pr_y2_povmed_v = sqrt(V[1,1])
scalar pr_y2_povmed_lci = pr_y2_povmed-invnormal(0.975)*pr_y2_povmed_v
scalar pr_y2_povmed_uci = pr_y2_povmed+invnormal(0.975)*pr_y2_povmed_v

disp "probability that y=2 is " pr_y2_povmed " and the 95% confidence interval is from "pr_y2_povmed_lci " to " pr_y2_povmed_uci
file write csvlog " 2 ," %20.5f (pr_y2_povmed) "," %20.5f (pr_y2_povmed_lci) "," %20.5f (pr_y2_povmed_uci) "," _n

** CT model (W= hood poverty) with covariates 

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_perpov_dw_z = x_f_release1 $x_covariates interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml init xb1:x_f_site_bal=0.1825 x_f_site_bos=0.0528 x_f_site_chi=0.1811 x_f_site_la=0.0674 f_c9010t_perpov_dw_z=-0.1610
ml maximize, difficult

*** counterfactual probabilities that Y = 0,1,2 (reported in column 2 of Table 6)

file write csvlog "y, P(Y=y|povX), lower CI, upper CI" _n

nlcom normal([cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}+ ${xgamma_morex}) )
matrix b = r(b)
matrix V = r(V)
scalar pr_x_y0_povmed = b[1,1]
scalar pr_x_y0_povmed_v = sqrt(V[1,1])
scalar pr_x_y0_povmed_lci = pr_x_y0_povmed-invnormal(0.975)*pr_x_y0_povmed_v
scalar pr_x_y0_povmed_uci = pr_x_y0_povmed+invnormal(0.975)*pr_x_y0_povmed_v

disp "probability that y=0 is " pr_x_y0_povmed " and the 95% confidence interval is from "pr_x_y0_povmed_lci " to " pr_x_y0_povmed_uci
file write csvlog " 0 ," %20.5f (pr_x_y0_povmed) "," %20.5f (pr_x_y0_povmed_lci) "," %20.5f (pr_x_y0_povmed_uci) "," _n

nlcom normal(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}+ ${xgamma_morex}) ) - /*
							*/normal([cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}+ ${xgamma_morex}) )
matrix b = r(b)
matrix V = r(V)
scalar pr_x_y1_povmed = b[1,1]
scalar pr_x_y1_povmed_v = sqrt(V[1,1])
scalar pr_x_y1_povmed_lci = pr_x_y1_povmed-invnormal(0.975)*pr_x_y1_povmed_v
scalar pr_x_y1_povmed_uci = pr_x_y1_povmed+invnormal(0.975)*pr_x_y1_povmed_v

disp "probability that y=1 is " pr_x_y1_povmed " and the 95% confidence interval is from "pr_x_y1_povmed_lci " to " pr_x_y1_povmed_uci
file write csvlog " 1 ," %20.5f (pr_x_y1_povmed) "," %20.5f (pr_x_y1_povmed_lci) "," %20.5f (pr_x_y1_povmed_uci) "," _n

nlcom 1-normal(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}+ ${xgamma_morex}) ) 
matrix b = r(b)
matrix V = r(V)
scalar pr_x_y2_povmed = b[1,1]
scalar pr_x_y2_povmed_v = sqrt(V[1,1])
scalar pr_x_y2_povmed_lci = pr_x_y2_povmed-invnormal(0.975)*pr_x_y2_povmed_v
scalar pr_x_y2_povmed_uci = pr_x_y2_povmed+invnormal(0.975)*pr_x_y2_povmed_v

disp "probability that y=2 is " pr_x_y2_povmed " and the 95% confidence interval is from "pr_x_y2_povmed_lci " to " pr_x_y2_povmed_uci
file write csvlog " 2 ," %20.5f (pr_x_y2_povmed) "," %20.5f (pr_x_y2_povmed_lci) "," %20.5f (pr_x_y2_povmed_uci) "," _n

file close csvlog

** TABLE 6b: columns 4 and 5 

file open csvlog using "${pathres}table6b_${date}.csv", write replace

** now re-set value of hood poverty at (median-SD) to estimate counterfactual probabilities reported in
** columns 4 and 5 of Table 6

scalar hood_pov_touse_nox = f_c9010t_perpov_dw_z_med_nox-f_c9010t_perpov_dw_z_sd_nox
scalar hood_pov_touse_x = f_c9010t_perpov_dw_z_med_x-f_c9010t_perpov_dw_z_sd_x

** CT probabilities 
** CT model (W=hood poverty) without covariates 

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs f_c9010t_perpov_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_perpov_dw_z = interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml init xb1:x_f_site_bal=0.1825 x_f_site_bos=0.0528 x_f_site_chi=0.1811 x_f_site_la=0.0674 f_c9010t_perpov_dw_z=-0.1610
ml maximize, difficult

*** counterfactual probabilities that Y = 0,1,2 (reported in column 4, Table 6)

file write csvlog "y, P(Y=y|pov), lower CI, upper CI" _n

nlcom normal([cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}))
matrix b = r(b)
matrix V = r(V)
scalar pr_y0_povmed = b[1,1]
scalar pr_y0_povmed_v = sqrt(V[1,1])
scalar pr_y0_povmed_lci = pr_y0_povmed-invnormal(0.975)*pr_y0_povmed_v
scalar pr_y0_povmed_uci = pr_y0_povmed+invnormal(0.975)*pr_y0_povmed_v

disp "probability that y=0 is " pr_y0_povmed " and the 95% confidence interval is from "pr_y0_povmed_lci " to " pr_y0_povmed_uci
file write csvlog " 0 ," %20.5f (pr_y0_povmed) "," %20.5f (pr_y0_povmed_lci) "," %20.5f (pr_y0_povmed_uci) "," _n

nlcom normal(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short})) - /*
							*/normal([cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}))
matrix b = r(b)
matrix V = r(V)
scalar pr_y1_povmed = b[1,1]
scalar pr_y1_povmed_v = sqrt(V[1,1])
scalar pr_y1_povmed_lci = pr_y1_povmed-invnormal(0.975)*pr_y1_povmed_v
scalar pr_y1_povmed_uci = pr_y1_povmed+invnormal(0.975)*pr_y1_povmed_v

disp "probability that y=1 is " pr_y1_povmed " and the 95% confidence interval is from "pr_y1_povmed_lci " to " pr_y1_povmed_uci
file write csvlog " 1 ," %20.5f (pr_y1_povmed) "," %20.5f (pr_y1_povmed_lci) "," %20.5f (pr_y1_povmed_uci) "," _n

nlcom 1-normal(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_nox)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short})) 
matrix b = r(b)
matrix V = r(V)
scalar pr_y2_povmed = b[1,1]
scalar pr_y2_povmed_v = sqrt(V[1,1])
scalar pr_y2_povmed_lci = pr_y2_povmed-invnormal(0.975)*pr_y2_povmed_v
scalar pr_y2_povmed_uci = pr_y2_povmed+invnormal(0.975)*pr_y2_povmed_v

disp "probability that y=2 is " pr_y2_povmed " and the 95% confidence interval is from "pr_y2_povmed_lci " to " pr_y2_povmed_uci
file write csvlog " 2 ," %20.5f (pr_y2_povmed) "," %20.5f (pr_y2_povmed_lci) "," %20.5f (pr_y2_povmed_uci) "," _n

** CT model (W= hood poverty) with covariates 

file write csvlog "y, P(Y=y|povX), lower CI, upper CI" _n

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_perpov_dw_z = x_f_release1 $x_covariates interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml init xb1:x_f_site_bal=0.1825 x_f_site_bos=0.0528 x_f_site_chi=0.1811 x_f_site_la=0.0674 f_c9010t_perpov_dw_z=-0.1610
ml maximize, difficult

*** counterfactual probabilities that Y = 0,1,2 (column 5, Table 6)

nlcom normal([cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}+ ${xgamma_morex}) )
matrix b = r(b)
matrix V = r(V)
scalar pr_x_y0_povmed = b[1,1]
scalar pr_x_y0_povmed_v = sqrt(V[1,1])
scalar pr_x_y0_povmed_lci = pr_x_y0_povmed-invnormal(0.975)*pr_x_y0_povmed_v
scalar pr_x_y0_povmed_uci = pr_x_y0_povmed+invnormal(0.975)*pr_x_y0_povmed_v

disp "probability that y=0 is " pr_x_y0_povmed " and the 95% confidence interval is from "pr_x_y0_povmed_lci " to " pr_x_y0_povmed_uci
file write csvlog " 0 ," %20.5f (pr_x_y0_povmed) "," %20.5f (pr_x_y0_povmed_lci) "," %20.5f (pr_x_y0_povmed_uci) "," _n

nlcom normal(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}+ ${xgamma_morex}) ) - /*
							*/normal([cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}+ ${xgamma_morex}) )
matrix b = r(b)
matrix V = r(V)
scalar pr_x_y1_povmed = b[1,1]
scalar pr_x_y1_povmed_v = sqrt(V[1,1])
scalar pr_x_y1_povmed_lci = pr_x_y1_povmed-invnormal(0.975)*pr_x_y1_povmed_v
scalar pr_x_y1_povmed_uci = pr_x_y1_povmed+invnormal(0.975)*pr_x_y1_povmed_v

disp "probability that y=1 is " pr_x_y1_povmed " and the 95% confidence interval is from "pr_x_y1_povmed_lci " to " pr_x_y1_povmed_uci
file write csvlog " 1 ," %20.5f (pr_x_y1_povmed) "," %20.5f (pr_x_y1_povmed_lci) "," %20.5f (pr_x_y1_povmed_uci) "," _n

nlcom 1-normal(exp([lndiff]_cons)+[cut1]_cons-(hood_pov_touse_x)*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short}+ ${xgamma_morex}) ) 
matrix b = r(b)
matrix V = r(V)
scalar pr_x_y2_povmed = b[1,1]
scalar pr_x_y2_povmed_v = sqrt(V[1,1])
scalar pr_x_y2_povmed_lci = pr_x_y2_povmed-invnormal(0.975)*pr_x_y2_povmed_v
scalar pr_x_y2_povmed_uci = pr_x_y2_povmed+invnormal(0.975)*pr_x_y2_povmed_v

disp "probability that y=2 is " pr_x_y2_povmed " and the 95% confidence interval is from "pr_x_y2_povmed_lci " to " pr_x_y2_povmed_uci
file write csvlog " 2 ," %20.5f (pr_x_y2_povmed) "," %20.5f (pr_x_y2_povmed_lci) "," %20.5f (pr_x_y2_povmed_uci) "," _n

file close csvlog
