
* generate variable medians for NY, used to get predicted probabilities and marginal effects later
 
sum f_c9010t_perpov_dw_z if (ra_site==5 & sample_nox==1), detail
scalar f_c9010t_perpov_dw_z_med_nox = r(p50)

sum f_c9010t_pminorty_dw_z if (ra_site==5 & sample_nox==1), detail
scalar f_c9010t_pminorty_dw_z_med_nox = r(p50)

sum f_c9010t_perpov_dw_z if (ra_site==5 & sample_x==1), detail
scalar f_c9010t_perpov_dw_z_med_x = r(p50)

sum f_c9010t_pminorty_dw_z if (ra_site==5 & sample_x==1), detail
scalar f_c9010t_pminorty_dw_z_med_x = r(p50)

foreach x in $x_covariates {
		
	summarize `x' if (ra_site==5 & sample_x==1), detail
	scalar `x'_med_x = r(p50)
	}
	
sum x_f_release1 if (ra_site==5 & sample_nox==1), detail
scalar x_f_release1_med_nox = r(p50)

sum x_f_release1 if (ra_site==5 & sample_x==1), detail
scalar x_f_release1_med_x = r(p50)

scalar x_rad_ad_41_45_med_x = 1
scalar x_f_hh_size3_med_x = 1

scalar x_rad_ad_ethrace_black_nh_med_x = 0 
scalar x_rad_ad_ethrace_hisp_med_x = 1 

** x_rad_ad_ethrace_black_nh_med = 0 
** x_rad_ad_ethrace_hisp_med = 1 
** tmp_female_med = 1
** x_rad_ad_le_35_med = 0 (9.7% of adults)
** x_rad_ad_36_40_med = 0 (20.7% of adults)
** x_rad_ad_41_45_med = 0 (force set to 1; 23.9% of adults)
** x_rad_ad_46_50_med = 0 (18.1% of adults)
** x_f_ad_nevmarr_med = 1
** x_f_ad_parentu18_med = 0
** x_f_ad_working_med = 0
** x_f_ad_edinsch_med = 0
** x_f_ad_edgradhs_med = 0
** x_f_ad_edged_med = 0
** x_f_hh_afdc_med = 1
** rad_svy_bl_totincm_2009d_med = 11950
** x_f_hh_car_med = 0
** x_f_hh_disabl_med = 0
** x_f_hh_noteens_med = 1
** x_f_hh_size2_med = 0 (22% of adults)
** x_f_hh_size3_med = 0 (force set to 1; 29% of adults)
** x_f_hh_size4_med = 0 (25% of adults)
** x_f_hh_victim_med = 1
** x_f_hood_unsafenit_med = 1
** x_f_hood_verydissat_med = 1
** x_f_hood_5y_med = 1
** x_f_hous_mov3tm_med = 0
** x_f_hood_nofamily_med = 1
** x_f_hood_nofriend_med = 0
** x_f_hood_chat_med = 1
** x_f_hood_nbrkid_med = 0
** x_f_hous_fndapt_med = 0
** x_f_hous_sec8bef_med = 0
** x_f_hous_movdrgs_med = 1
** x_f_hous_movschl_med = 1
** cov_hous_movapt_med = 0
** cov_hous_movjob_med = 0
** x_f_release1_med = 0

* expressions to be used later for predicted probabilities and marginal effects

global x_fix_short "x_f_site_balt=0 x_f_site_bos=0 x_f_site_chi=0 x_f_site_la=0"
global x_fix_long  "x_rad_ad_ethrace_black_nh=0 x_rad_ad_ethrace_hisp=1 tmp_female=1 x_rad_ad_le_35=0 x_rad_ad_36_40=0 x_rad_ad_41_45=1 x_rad_ad_46_50=0 x_f_ad_nevmarr=1 x_f_ad_parentu18=0 x_f_ad_working=0 x_f_ad_edinsch=0 x_f_ad_edgradhs=0 x_f_ad_edged=0 x_f_hh_afdc=1 rad_svy_bl_totincm_2009d=11950 x_f_hh_car=0 x_f_hh_disabl=0 x_f_hh_noteens=1 x_f_hh_size2=0 x_f_hh_size3=1 x_f_hh_size4=0 x_f_hh_victim=1 x_f_hood_unsafenit=1 x_f_hood_verydissat=1 x_f_hood_5y=1 x_f_hous_mov3tm=0 x_f_hood_nofamily=1 x_f_hood_nofriend=0 x_f_hood_chat=1 x_f_hood_nbrkid=0 x_f_hous_fndapt=0 x_f_hous_sec8bef=0 x_f_hous_movdrgs=1 x_f_hous_movschl=1 cov_hous_movapt=0 cov_hous_movjob=0 x_f_release1=0 x_f_site_balt=0 x_f_site_bos=0 x_f_site_chi=0 x_f_site_la=0"

global xgamma_short "(0) "
global xgamma_morex "x_f_release1_med_x*[xb1]x_f_release1+tmp_female_med_x*[xb1]tmp_female+x_rad_ad_le_35_med_x*[xb1]x_rad_ad_le_35+x_rad_ad_36_40_med_x*[xb1]x_rad_ad_36_40+x_rad_ad_41_45_med_x*[xb1]x_rad_ad_41_45+x_rad_ad_46_50_med_x*[xb1]x_rad_ad_46_50+x_rad_ad_ethrace_black_nh_med_x*[xb1]x_rad_ad_ethrace_black_nh+x_rad_ad_ethrace_hisp_med_x*[xb1]x_rad_ad_ethrace_hisp+x_f_ad_nevmarr_med_x*[xb1]x_f_ad_nevmarr+x_f_ad_parentu18_med_x*[xb1]x_f_ad_parentu18+x_f_ad_working_med_x*[xb1]x_f_ad_working+x_f_ad_edinsch_med_x*[xb1]x_f_ad_edinsch+x_f_ad_edgradhs_med_x*[xb1]x_f_ad_edgradhs+x_f_ad_edged_med_x*[xb1]x_f_ad_edged+x_f_hh_afdc_med_x*[xb1]x_f_hh_afdc+rad_svy_bl_totincm_2009d_med_x*[xb1]rad_svy_bl_totincm_2009d+x_f_hh_car_med_x*[xb1]x_f_hh_car+x_f_hh_disabl_med_x*[xb1]x_f_hh_disabl+x_f_hh_noteens_med_x*[xb1]x_f_hh_noteens+x_f_hh_size2_med_x*[xb1]x_f_hh_size2+x_f_hh_size3_med_x*[xb1]x_f_hh_size3+x_f_hh_size4_med_x*[xb1]x_f_hh_size4+x_f_hh_victim_med_x*[xb1]x_f_hh_victim+x_f_hood_unsafenit_med_x*[xb1]x_f_hood_unsafenit+x_f_hood_verydissat_med_x*[xb1]x_f_hood_verydissat+x_f_hood_5y_med_x*[xb1]x_f_hood_5y+x_f_hous_mov3tm_med_x*[xb1]x_f_hous_mov3tm+x_f_hood_nofamily_med_x*[xb1]x_f_hood_nofamily+x_f_hood_nofriend_med_x*[xb1]x_f_hood_nofriend+x_f_hood_chat_med_x*[xb1]x_f_hood_chat+x_f_hood_nbrkid_med_x*[xb1]x_f_hood_nbrkid+x_f_hous_fndapt_med_x*[xb1]x_f_hous_fndapt+x_f_hous_sec8bef_med_x*[xb1]x_f_hous_sec8bef+x_f_hous_movdrgs_med_x*[xb1]x_f_hous_movdrgs+x_f_hous_movschl_med_x*[xb1]x_f_hous_movschl+cov_hous_movapt_med_x*[xb1]cov_hous_movapt+cov_hous_movjob_med_x*[xb1]cov_hous_movjob"

** ordered probit probabilities and marginal effects

** FIGURE 3 

** W = hood poverty, with covariates

** CT model 

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_perpov_dw_z = x_f_release1 $x_covariates interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml init xb1:x_f_site_bal=0.1825 x_f_site_bos=0.0528 x_f_site_chi=0.1811 x_f_site_la=0.0674 f_c9010t_perpov_dw_z=-0.1610
ml maximize, difficult
est store pov_t_x

forvalues y=0/2{
	forvalues i = 1/11{
		scalar poverty_`y'_`i' = `i' - 7
		if `y'==0{
		nlcom normal([cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short} + ${xgamma_morex}) )
		}
		else if `y'==1{
		nlcom normal(exp([lndiff]_cons)+[cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short} + ${xgamma_morex}) ) - /*				
				*/normal([cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short} + ${xgamma_morex}) )
		}
		else if `y'==2{
		nlcom 1-normal(exp([lndiff]_cons)+[cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short} + ${xgamma_morex}) ) 
		}
		matrix b = r(b)
		matrix V = r(V)
		scalar py`y'_`i'_med_ny = b[1,1]
		scalar v`y'_`i' = sqrt(V[1,1])
		scalar lci`y'_`i'_med_ny = py`y'_`i'_med_ny-invnormal(0.975)*v`y'_`i'
		scalar uci`y'_`i'_med_ny = py`y'_`i'_med_ny+invnormal(0.975)*v`y'_`i'
	}

	matrix py`y'_x_pov_med_ny = (scalar(poverty_`y'_1),scalar(py`y'_1_med_ny),scalar(lci`y'_1_med_ny),scalar(uci`y'_1_med_ny)\/*
				*/scalar(poverty_`y'_2),scalar(py`y'_2_med_ny),scalar(lci`y'_2_med_ny),scalar(uci`y'_2_med_ny)\/*
				*/scalar(poverty_`y'_3),scalar(py`y'_3_med_ny),scalar(lci`y'_3_med_ny),scalar(uci`y'_3_med_ny)\/*
				*/scalar(poverty_`y'_4),scalar(py`y'_4_med_ny),scalar(lci`y'_4_med_ny),scalar(uci`y'_4_med_ny)\/*
				*/scalar(poverty_`y'_5),scalar(py`y'_5_med_ny),scalar(lci`y'_5_med_ny),scalar(uci`y'_5_med_ny)\/*
				*/scalar(poverty_`y'_6),scalar(py`y'_6_med_ny),scalar(lci`y'_6_med_ny),scalar(uci`y'_6_med_ny)\/*
				*/scalar(poverty_`y'_7),scalar(py`y'_7_med_ny),scalar(lci`y'_7_med_ny),scalar(uci`y'_7_med_ny)\/*
				*/scalar(poverty_`y'_8),scalar(py`y'_8_med_ny),scalar(lci`y'_8_med_ny),scalar(uci`y'_8_med_ny)\/*
				*/scalar(poverty_`y'_9),scalar(py`y'_9_med_ny),scalar(lci`y'_9_med_ny),scalar(uci`y'_9_med_ny)\/*
				*/scalar(poverty_`y'_10),scalar(py`y'_10_med_ny),scalar(lci`y'_10_med_ny),scalar(uci`y'_10_med_ny)\/*
				*/scalar(poverty_`y'_11),scalar(py`y'_11_med_ny),scalar(lci`y'_11_med_ny),scalar(uci`y'_11_med_ny))

	matrix colnames py`y'_x_pov_med_ny = povx`y'_med_ny pyx`y'_med_ny lcix`y'_med_ny ucix`y'_med_ny
	svmat py`y'_x_pov_med_ny, names(col)
	
	matrix list py`y'_x_pov_med_ny
}

** ordered probit model 

ml model lf myoprobit_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z, nocons) (cut1:) (cut2:) [pw=$wt]
ml maximize
est store pov_op_x

forvalues y = 0/2 {
	if `y'==0{
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=(-6(1)4)) expression(normal(predict(equation(cut1))-predict(equation(xb1))))	
	}
	else if `y'==1{
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=(-6(1)4)) expression((normal(predict(equation(cut2))-predict(equation(xb1))))-(normal(predict(equation(cut1))-predict(equation(xb1)))))
	}
	else if `y'==2{
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=(-6(1)4)) expression(1-(normal(predict(equation(cut2))-predict(equation(xb1)))))
	}
	
	marginsplot, ylabel(-0.4(0.4)1.2, grid) title(" ") xtitle(" ") ytitle(" ") xlabel(-6(2)4) recast(line) recastci(rarea) plotopts(lcolor(black) lwidth(medthick)) ciopts(color(black%50)) ///
	addplot(rarea lcix`y'_med_ny ucix`y'_med_ny povx`y'_med_ny, ylabel(-0.4(0.4)1.2, grid) color(black%20)||line pyx`y'_med_ny povx`y'_med_ny, ylabel(-0.4(0.4)1.2, grid) ///
	lcolor(gs5) lpattern(dash) lwidth(medthick)) legend(off)
	graph export ${pathres}prob_y`y'_x_pov_med_ny_${date}.pdf, replace
	}
	
** W = hood minority, with covariates

** CT model 

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_pminorty_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_pminorty_dw_z = x_f_release1 $x_covariates interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml maximize
est store min_t_x

forvalues y = 0/2{
	forvalues i = 1/11{
		scalar minority_`y'_`i' = `i' - 7
		if `y'==0{
		nlcom normal([cut1]_cons-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z - (${xgamma_short} + ${xgamma_morex}) )
		}
		else if `y'==1{
		nlcom normal(exp([lndiff]_cons)+[cut1]_cons-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z - (${xgamma_short} + ${xgamma_morex}) ) - /*				
				*/normal([cut1]_cons-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z - (${xgamma_short} + ${xgamma_morex}) )
		}
		else if `y'==2{
		nlcom 1-normal(exp([lndiff]_cons)+[cut1]_cons-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z - (${xgamma_short} + ${xgamma_morex}) ) 
		}
		
		matrix b = r(b)
		matrix V = r(V)
		scalar py`y'_`i'_med_ny = b[1,1]
		scalar v`y'_`i' = sqrt(V[1,1])
		scalar lci`y'_`i'_med_ny = py`y'_`i'_med_ny-invnormal(0.975)*v`y'_`i'
		scalar uci`y'_`i'_med_ny = py`y'_`i'_med_ny+invnormal(0.975)*v`y'_`i'
	}

	matrix py`y'_x_min_med_ny = (scalar(minority_`y'_1),scalar(py`y'_1_med_ny),scalar(lci`y'_1_med_ny),scalar(uci`y'_1_med_ny)\/*
				*/scalar(minority_`y'_2),scalar(py`y'_2_med_ny),scalar(lci`y'_2_med_ny),scalar(uci`y'_2_med_ny)\/*
				*/scalar(minority_`y'_3),scalar(py`y'_3_med_ny),scalar(lci`y'_3_med_ny),scalar(uci`y'_3_med_ny)\/*
				*/scalar(minority_`y'_4),scalar(py`y'_4_med_ny),scalar(lci`y'_4_med_ny),scalar(uci`y'_4_med_ny)\/*
				*/scalar(minority_`y'_5),scalar(py`y'_5_med_ny),scalar(lci`y'_5_med_ny),scalar(uci`y'_5_med_ny)\/*
				*/scalar(minority_`y'_6),scalar(py`y'_6_med_ny),scalar(lci`y'_6_med_ny),scalar(uci`y'_6_med_ny)\/*
				*/scalar(minority_`y'_7),scalar(py`y'_7_med_ny),scalar(lci`y'_7_med_ny),scalar(uci`y'_7_med_ny)\/*
				*/scalar(minority_`y'_8),scalar(py`y'_8_med_ny),scalar(lci`y'_8_med_ny),scalar(uci`y'_8_med_ny)\/*
				*/scalar(minority_`y'_9),scalar(py`y'_9_med_ny),scalar(lci`y'_9_med_ny),scalar(uci`y'_9_med_ny)\/*
				*/scalar(minority_`y'_10),scalar(py`y'_10_med_ny),scalar(lci`y'_10_med_ny),scalar(uci`y'_10_med_ny)\/*
				*/scalar(minority_`y'_11),scalar(py`y'_11_med_ny),scalar(lci`y'_11_med_ny),scalar(uci`y'_11_med_ny))

	matrix colnames py`y'_x_min_med_ny = minx`y'_med_ny minpyx`y'_med_ny minlcix`y'_med_ny minucix`y'_med_ny
	svmat py`y'_x_min_med_ny, names(col)
	
	matrix list py`y'_x_min_med_ny
}

** ordered probit model 

ml model lf myoprobit_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_pminorty_dw_z, nocons) (cut1:) (cut2:) [pw=$wt]
ml maximize
est store min_op_x

forvalues y = 0/2 {
	if `y'==0{
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=(-6(1)4)) expression(normal(predict(equation(cut1))-predict(equation(xb1))))
	}
	else if `y'==1{
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=(-6(1)4)) expression((normal(predict(equation(cut2))-predict(equation(xb1))))-(normal(predict(equation(cut1))-predict(equation(xb1)))))
	}
	else if `y'==2{
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=(-6(1)4)) expression(1-(normal(predict(equation(cut2))-predict(equation(xb1)))))
	}
	
	marginsplot, ylabel(-0.4(0.4)1.2, grid) title(" ") xtitle(" ") ytitle(" ") xlabel(-6(2)4) recast(line) recastci(rarea) plotopts(lcolor(black) lwidth(medthick)) ciopts(color(black%50)) ///
	addplot(rarea minlcix`y'_med_ny minucix`y'_med_ny minx`y'_med_ny, ylabel(-0.4(0.4)1.2, grid) color(black%20)||line minpyx`y'_med_ny minx`y'_med_ny, ylabel(-0.4(0.4)1.2, grid) ///
	lcolor(gs5) lpattern(dash) lwidth(medthick)) legend(off)
	graph export ${pathres}prob_y`y'_x_min_med_ny_${date}.pdf, replace
	}

** FIGURE 4 
	
** W = (hood poverty, hood minority), with covariates

** CT model 

ml model lf myoprobit3_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z, nocons) (cut1:) (lndiff:) /*
						*/(xb2:f_c9010t_perpov_dw_z = x_f_release1 $x_covariates interaction*) (lnvar_v1:) (rho_uv1:) /*
						*/(xb3:f_c9010t_pminorty_dw_z = x_f_release1 $x_covariates interaction*) (lnvar_v2:) (rho_uv2:) (cov_v1v2:) [pw=$wt], /*
						*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
						*/diparm(lnvar_v1, exp label(var_v1)) /*
						*/diparm(lnvar_v2, exp label(var_v2)) 
ml init xb1:x_f_site_bal=0.3039 x_f_site_bos=0.3386 x_f_site_chi=0.1527 x_f_site_la=0.0639 f_c9010t_perpov_dw_z=-0.2998 f_c9010t_pminorty_dw_z=0.3151 cut1:_cons=-0.8165
ml maximize, difficult
est store povmin_t_x

forvalues y = 0/2{
	forvalues i = 1/11{
		scalar poverty_`y'_`i' = `i' - 7
		if `y'==0{
		nlcom normal([cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_x*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))
		matrix b = r(b)
		matrix V = r(V)
		scalar povmin1_x_py`y'_`i'_med_ny = b[1,1]
		scalar povmin1_x_v`y'_`i' = sqrt(V[1,1])
		scalar povmin1_x_lci`y'_`i'_med_ny = povmin1_x_py`y'_`i'_med_ny-invnormal(0.975)*povmin1_x_v`y'_`i'
		scalar povmin1_x_uci`y'_`i'_med_ny = povmin1_x_py`y'_`i'_med_ny+invnormal(0.975)*povmin1_x_v`y'_`i'
		}
		else if `y'==1{
		nlcom normal(exp([lndiff]_cons)+[cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_x*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex})) - /*
				*/normal([cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_x*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))
		matrix b = r(b)
		matrix V = r(V)
		scalar povmin1_x_py`y'_`i'_med_ny = b[1,1]
		scalar povmin1_x_v`y'_`i' = sqrt(V[1,1])
		scalar povmin1_x_lci`y'_`i'_med_ny = povmin1_x_py`y'_`i'_med_ny-invnormal(0.975)*povmin1_x_v`y'_`i'
		scalar povmin1_x_uci`y'_`i'_med_ny = povmin1_x_py`y'_`i'_med_ny+invnormal(0.975)*povmin1_x_v`y'_`i'
		}
		else if `y'==2{
		nlcom 1-normal(exp([lndiff]_cons)+[cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_x*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))
		matrix b = r(b)
		matrix V = r(V)
		scalar povmin1_x_py`y'_`i'_med_ny = b[1,1]
		scalar povmin1_x_v`y'_`i' = sqrt(V[1,1])
		scalar povmin1_x_lci`y'_`i'_med_ny = povmin1_x_py`y'_`i'_med_ny-invnormal(0.975)*povmin1_x_v`y'_`i'
		scalar povmin1_x_uci`y'_`i'_med_ny = povmin1_x_py`y'_`i'_med_ny+invnormal(0.975)*povmin1_x_v`y'_`i'
		}
	}

	matrix py`y'_x_povmin1_med_ny = (scalar(poverty_`y'_1),scalar(povmin1_x_py`y'_1_med_ny),scalar(povmin1_x_lci`y'_1_med_ny),scalar(povmin1_x_uci`y'_1_med_ny)\/*
				*/scalar(poverty_`y'_2),scalar(povmin1_x_py`y'_2_med_ny),scalar(povmin1_x_lci`y'_2_med_ny),scalar(povmin1_x_uci`y'_2_med_ny)\/*
				*/scalar(poverty_`y'_3),scalar(povmin1_x_py`y'_3_med_ny),scalar(povmin1_x_lci`y'_3_med_ny),scalar(povmin1_x_uci`y'_3_med_ny)\/*
				*/scalar(poverty_`y'_4),scalar(povmin1_x_py`y'_4_med_ny),scalar(povmin1_x_lci`y'_4_med_ny),scalar(povmin1_x_uci`y'_4_med_ny)\/*
				*/scalar(poverty_`y'_5),scalar(povmin1_x_py`y'_5_med_ny),scalar(povmin1_x_lci`y'_5_med_ny),scalar(povmin1_x_uci`y'_5_med_ny)\/*
				*/scalar(poverty_`y'_6),scalar(povmin1_x_py`y'_6_med_ny),scalar(povmin1_x_lci`y'_6_med_ny),scalar(povmin1_x_uci`y'_6_med_ny)\/*
				*/scalar(poverty_`y'_7),scalar(povmin1_x_py`y'_7_med_ny),scalar(povmin1_x_lci`y'_7_med_ny),scalar(povmin1_x_uci`y'_7_med_ny)\/*
				*/scalar(poverty_`y'_8),scalar(povmin1_x_py`y'_8_med_ny),scalar(povmin1_x_lci`y'_8_med_ny),scalar(povmin1_x_uci`y'_8_med_ny)\/*
				*/scalar(poverty_`y'_9),scalar(povmin1_x_py`y'_9_med_ny),scalar(povmin1_x_lci`y'_9_med_ny),scalar(povmin1_x_uci`y'_9_med_ny)\/*
				*/scalar(poverty_`y'_10),scalar(povmin1_x_py`y'_10_med_ny),scalar(povmin1_x_lci`y'_10_med_ny),scalar(povmin1_x_uci`y'_10_med_ny)\/*
				*/scalar(poverty_`y'_11),scalar(povmin1_x_py`y'_11_med_ny),scalar(povmin1_x_lci`y'_11_med_ny),scalar(povmin1_x_uci`y'_11_med_ny))

	matrix colnames py`y'_x_povmin1_med_ny = povmin1_x_pov`y'_med_ny povmin1_x_py`y'_med_ny povmin1_x_lci`y'_med_ny povmin1_x_uci`y'_med_ny
	svmat py`y'_x_povmin1_med_ny, names(col)

	matrix list py`y'_x_povmin1_med_ny
}

forvalues y = 0/2{
	forvalues i = 1/11{
		scalar minority_`y'_`i' = `i' - 7
		if `y'==0{
		nlcom normal([cut1]_cons-f_c9010t_perpov_dw_z_med_x*[xb1]f_c9010t_perpov_dw_z-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))
		matrix b = r(b)
		matrix V = r(V)
		scalar povmin2_x_py`y'_`i'_med_ny = b[1,1]
		scalar povmin2_x_v`y'_`i' = sqrt(V[1,1])
		scalar povmin2_x_lci`y'_`i'_med_ny = povmin2_x_py`y'_`i'_med_ny-invnormal(0.975)*povmin2_x_v`y'_`i'
		scalar povmin2_x_uci`y'_`i'_med_ny = povmin2_x_py`y'_`i'_med_ny+invnormal(0.975)*povmin2_x_v`y'_`i'
		}
		else if `y'==1{
		nlcom normal(exp([lndiff]_cons)+[cut1]_cons-f_c9010t_perpov_dw_z_med_x*[xb1]f_c9010t_perpov_dw_z-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex})) - /*
				*/normal([cut1]_cons-f_c9010t_perpov_dw_z_med_x*[xb1]f_c9010t_perpov_dw_z-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))
		matrix b = r(b)
		matrix V = r(V)
		scalar povmin2_x_py`y'_`i'_med_ny = b[1,1]
		scalar povmin2_x_v`y'_`i' = sqrt(V[1,1])
		scalar povmin2_x_lci`y'_`i'_med_ny = povmin2_x_py`y'_`i'_med_ny-invnormal(0.975)*povmin2_x_v`y'_`i'
		scalar povmin2_x_uci`y'_`i'_med_ny = povmin2_x_py`y'_`i'_med_ny+invnormal(0.975)*povmin2_x_v`y'_`i'
		}
		else if `y'==2{
		nlcom 1-normal(exp([lndiff]_cons)+[cut1]_cons-f_c9010t_perpov_dw_z_med_x*[xb1]f_c9010t_perpov_dw_z-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))
		matrix b = r(b)
		matrix V = r(V)
		scalar povmin2_x_py`y'_`i'_med_ny = b[1,1]
		scalar povmin2_x_v`y'_`i' = sqrt(V[1,1])
		scalar povmin2_x_lci`y'_`i'_med_ny = povmin2_x_py`y'_`i'_med_ny-invnormal(0.975)*povmin2_x_v`y'_`i'
		scalar povmin2_x_uci`y'_`i'_med_ny = povmin2_x_py`y'_`i'_med_ny+invnormal(0.975)*povmin2_x_v`y'_`i'
		}
	}

	matrix py`y'_x_povmin2_med_ny = (scalar(minority_`y'_1),scalar(povmin2_x_py`y'_1_med_ny),scalar(povmin2_x_lci`y'_1_med_ny),scalar(povmin2_x_uci`y'_1_med_ny)\/*
				*/scalar(minority_`y'_2),scalar(povmin2_x_py`y'_2_med_ny),scalar(povmin2_x_lci`y'_2_med_ny),scalar(povmin2_x_uci`y'_2_med_ny)\/*
				*/scalar(minority_`y'_3),scalar(povmin2_x_py`y'_3_med_ny),scalar(povmin2_x_lci`y'_3_med_ny),scalar(povmin2_x_uci`y'_3_med_ny)\/*
				*/scalar(minority_`y'_4),scalar(povmin2_x_py`y'_4_med_ny),scalar(povmin2_x_lci`y'_4_med_ny),scalar(povmin2_x_uci`y'_4_med_ny)\/*
				*/scalar(minority_`y'_5),scalar(povmin2_x_py`y'_5_med_ny),scalar(povmin2_x_lci`y'_5_med_ny),scalar(povmin2_x_uci`y'_5_med_ny)\/*
				*/scalar(minority_`y'_6),scalar(povmin2_x_py`y'_6_med_ny),scalar(povmin2_x_lci`y'_6_med_ny),scalar(povmin2_x_uci`y'_6_med_ny)\/*
				*/scalar(minority_`y'_7),scalar(povmin2_x_py`y'_7_med_ny),scalar(povmin2_x_lci`y'_7_med_ny),scalar(povmin2_x_uci`y'_7_med_ny)\/*
				*/scalar(minority_`y'_8),scalar(povmin2_x_py`y'_8_med_ny),scalar(povmin2_x_lci`y'_8_med_ny),scalar(povmin2_x_uci`y'_8_med_ny)\/*
				*/scalar(minority_`y'_9),scalar(povmin2_x_py`y'_9_med_ny),scalar(povmin2_x_lci`y'_9_med_ny),scalar(povmin2_x_uci`y'_9_med_ny)\/*
				*/scalar(minority_`y'_10),scalar(povmin2_x_py`y'_10_med_ny),scalar(povmin2_x_lci`y'_10_med_ny),scalar(povmin2_x_uci`y'_10_med_ny)\/*
				*/scalar(minority_`y'_11),scalar(povmin2_x_py`y'_11_med_ny),scalar(povmin2_x_lci`y'_11_med_ny),scalar(povmin2_x_uci`y'_11_med_ny))

	matrix colnames py`y'_x_povmin2_med_ny = povmin2_x_min`y'_med_ny povmin2_x_py`y'_med_ny povmin2_x_lci`y'_med_ny povmin2_x_uci`y'_med_ny
	svmat py`y'_x_povmin2_med_ny, names(col)
	
	matrix list py`y'_x_povmin2_med_ny
}

** ordered probit model 
	
ml model lf myoprobit_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z, nocons) (cut1:) (cut2:) [pw=$wt]
ml maximize
est store povmin_op_x

forvalues y = 0/2 {
	if `y'==0{
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=.4965963 f_c9010t_perpov_dw_z=(-6(1)4)) expression(normal(predict(equation(cut1))-predict(equation(xb1))))
	}
	else if `y'==1{
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=.4965963 f_c9010t_perpov_dw_z=(-6(1)4)) expression((normal(predict(equation(cut2))-predict(equation(xb1))))-(normal(predict(equation(cut1))-predict(equation(xb1)))))
	}
	else if `y'==2{
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=.4965963 f_c9010t_perpov_dw_z=(-6(1)4)) expression(1-(normal(predict(equation(cut2))-predict(equation(xb1)))))
	}
	
	marginsplot, ylabel(-0.6(0.3)1.5, grid) title(" ") xtitle(" ") ytitle(" ") xlabel(-6(2)4) recast(line) recastci(rarea) plotopts(lcolor(black) lwidth(medthick)) ciopts(color(black%50)) ///
	addplot(rarea povmin1_x_lci`y'_med_ny povmin1_x_uci`y'_med_ny povmin1_x_pov`y'_med_ny, ylabel(-0.6(0.3)1.5, grid) color(black%20)||line povmin1_x_py`y'_med_ny povmin1_x_pov`y'_med_ny, ///
	ylabel(-0.6(0.3)1.5, grid) lcolor(gs5) lpattern(dash) lwidth(medthick)) legend(off)
	graph export ${pathres}prob_y`y'_x_povmin1_med_ny_${date}.pdf, replace
	}
	
forvalues y = 0/2 {
	if `y'==0{
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=-.182967 f_c9010t_pminorty_dw_z=(-6(1)4)) expression(normal(predict(equation(cut1))-predict(equation(xb1))))
	}
	else if `y'==1{
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=-.182967 f_c9010t_pminorty_dw_z=(-6(1)4)) expression((normal(predict(equation(cut2))-predict(equation(xb1))))-(normal(predict(equation(cut1))-predict(equation(xb1)))))
	}
	else if `y'==2{
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=-.182967 f_c9010t_pminorty_dw_z=(-6(1)4)) expression(1-(normal(predict(equation(cut2))-predict(equation(xb1)))))
	}
	
	marginsplot, ylabel(-0.6(0.3)1.5, grid) title(" ") xtitle(" ") ytitle(" ") xlabel(-6(2)4) recast(line) recastci(rarea) plotopts(lcolor(black) lwidth(medthick)) ciopts(color(black%50)) ///
	addplot(rarea povmin2_x_lci`y'_med_ny povmin2_x_uci`y'_med_ny povmin2_x_min`y'_med_ny, ylabel(-0.6(0.3)1.5, grid) color(black%20)||line povmin2_x_py`y'_med_ny povmin2_x_min`y'_med_ny, ///
	ylabel(-0.6(0.3)1.5, grid) lcolor(gs5) lpattern(dash) lwidth(medthick)) legend(off)
	graph export ${pathres}prob_y`y'_x_povmin2_med_ny_${date}.pdf, replace
	}

** FIGURE 5 

** CT model 

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_perpov_dw_z = x_f_release1 $x_covariates interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml init xb1:x_f_site_bal=0.1825 x_f_site_bos=0.0528 x_f_site_chi=0.1811 x_f_site_la=0.0674 f_c9010t_perpov_dw_z=-0.1610
ml maximize, difficult
est store pov_t_x

forvalues y = 0/2 {
	forvalues i = 1/11{
		scalar poverty_`y'_`i' = `i' - 7
		if `y'==0{
		nlcom -[xb1]f_c9010t_perpov_dw_z * normalden([cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short} + ${xgamma_morex}) )
		}
		else if `y'==1{
		nlcom [xb1]f_c9010t_perpov_dw_z * (normalden([cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short} + ${xgamma_morex}) )- /*
				*/normalden(exp([lndiff]_cons)+[cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short} + ${xgamma_morex}) ))
		}
		else if `y'==2{
		nlcom [xb1]f_c9010t_perpov_dw_z * normalden(exp([lndiff]_cons)+[cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z - (${xgamma_short} + ${xgamma_morex}) )
		}
		matrix b = r(b)
		matrix V = r(V)
		scalar me`y'_`i'_med_ny = b[1,1]
		scalar v`y'_`i' = sqrt(V[1,1])
		scalar lci`y'_`i'_med_ny = me`y'_`i'_med_ny-invnormal(0.975)*v`y'_`i'
		scalar uci`y'_`i'_med_ny = me`y'_`i'_med_ny+invnormal(0.975)*v`y'_`i'
	}

	matrix me`y'_x_pov_med_ny = (scalar(poverty_`y'_1),scalar(me`y'_1_med_ny),scalar(lci`y'_1_med_ny),scalar(uci`y'_1_med_ny)\/*
				*/scalar(poverty_`y'_2),scalar(me`y'_2_med_ny),scalar(lci`y'_2_med_ny),scalar(uci`y'_2_med_ny)\/*
				*/scalar(poverty_`y'_3),scalar(me`y'_3_med_ny),scalar(lci`y'_3_med_ny),scalar(uci`y'_3_med_ny)\/*
				*/scalar(poverty_`y'_4),scalar(me`y'_4_med_ny),scalar(lci`y'_4_med_ny),scalar(uci`y'_4_med_ny)\/*
				*/scalar(poverty_`y'_5),scalar(me`y'_5_med_ny),scalar(lci`y'_5_med_ny),scalar(uci`y'_5_med_ny)\/*
				*/scalar(poverty_`y'_6),scalar(me`y'_6_med_ny),scalar(lci`y'_6_med_ny),scalar(uci`y'_6_med_ny)\/*
				*/scalar(poverty_`y'_7),scalar(me`y'_7_med_ny),scalar(lci`y'_7_med_ny),scalar(uci`y'_7_med_ny)\/*
				*/scalar(poverty_`y'_8),scalar(me`y'_8_med_ny),scalar(lci`y'_8_med_ny),scalar(uci`y'_8_med_ny)\/*
				*/scalar(poverty_`y'_9),scalar(me`y'_9_med_ny),scalar(lci`y'_9_med_ny),scalar(uci`y'_9_med_ny)\/*
				*/scalar(poverty_`y'_10),scalar(me`y'_10_med_ny),scalar(lci`y'_10_med_ny),scalar(uci`y'_10_med_ny)\/*
				*/scalar(poverty_`y'_11),scalar(me`y'_11_med_ny),scalar(lci`y'_11_med_ny),scalar(uci`y'_11_med_ny))

	matrix colnames me`y'_x_pov_med_ny = mepovx`y'_med_ny mex`y'_med_ny melcix`y'_med_ny meucix`y'_med_ny
	svmat me`y'_x_pov_med_ny, names(col)
	
	matrix list me`y'_x_pov_med_ny
}

** ordered probit model 

ml model lf myoprobit_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z, nocons) (cut1:) (cut2:) [pw=$wt]
ml maximize
est store pov_op_x

forvalues y = 0/2 {
	if `y'==0{
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=(-6(1)4)) dydx(f_c9010t_perpov_dw_z) expression(normal(predict(equation(cut1))-predict(equation(xb1))))
	}
	else if `y'==1{
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=(-6(1)4)) dydx(f_c9010t_perpov_dw_z) expression((normal(predict(equation(cut2))-predict(equation(xb1))))-(normal(predict(equation(cut1))-predict(equation(xb1)))))
	}
	else if `y'==2{
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=(-6(1)4)) dydx(f_c9010t_perpov_dw_z) expression(1-(normal(predict(equation(cut2))-predict(equation(xb1)))))
	}
	
	marginsplot, ylabel(-0.2(0.1)0.2, grid) title(" ") xtitle(" ") ytitle(" ") xlabel(-6(2)4) recast(line) recastci(rarea) plotopts(lcolor(black) lwidth(medthick)) ciopts(color(black%50)) ///
	addplot(rarea melcix`y'_med_ny meucix`y'_med_ny mepovx`y'_med_ny, ylabel(-0.2(0.1)0.2, grid) color(black%20)||line mex`y'_med_ny mepovx`y'_med_ny, ylabel(-0.2(0.1)0.2, grid) ///
	lcolor(gs5) lpattern(dash) lwidth(medthick)) legend(off)
	graph export ${pathres}me_y`y'_x_pov_med_ny_${date}.pdf, replace
	}

** W = hood minority, with covariates

** CT model 

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_pminorty_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_pminorty_dw_z = x_f_release1 $x_covariates interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml maximize
est store min_t_x

forvalues y = 0/2{
	forvalues i = 1/11{
		scalar minority_`y'_`i' = `i' - 7
		if `y'==0{
		nlcom -[xb1]f_c9010t_pminorty_dw_z * normalden([cut1]_cons-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z - (${xgamma_short} + ${xgamma_morex}) )
		}
		else if `y'==1{
		nlcom [xb1]f_c9010t_pminorty_dw_z * (normalden([cut1]_cons-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z - (${xgamma_short} + ${xgamma_morex}) )- /*
				*/normalden(exp([lndiff]_cons)+[cut1]_cons-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z - (${xgamma_short} + ${xgamma_morex}) ))
		}
		else if `y'==2{
		nlcom [xb1]f_c9010t_pminorty_dw_z * normalden(exp([lndiff]_cons)+[cut1]_cons-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z - (${xgamma_short} + ${xgamma_morex}) )
		}
		matrix b = r(b)
		matrix V = r(V)
		scalar me`y'_`i'_med_ny = b[1,1]
		scalar v`y'_`i' = sqrt(V[1,1])
		scalar lci`y'_`i'_med_ny = me`y'_`i'_med_ny-invnormal(0.975)*v`y'_`i'
		scalar uci`y'_`i'_med_ny = me`y'_`i'_med_ny+invnormal(0.975)*v`y'_`i'
	}

	matrix me`y'_x_min_med_ny = (scalar(minority_`y'_1),scalar(me`y'_1_med_ny),scalar(lci`y'_1_med_ny),scalar(uci`y'_1_med_ny)\/*
				*/scalar(minority_`y'_2),scalar(me`y'_2_med_ny),scalar(lci`y'_2_med_ny),scalar(uci`y'_2_med_ny)\/*
				*/scalar(minority_`y'_3),scalar(me`y'_3_med_ny),scalar(lci`y'_3_med_ny),scalar(uci`y'_3_med_ny)\/*
				*/scalar(minority_`y'_4),scalar(me`y'_4_med_ny),scalar(lci`y'_4_med_ny),scalar(uci`y'_4_med_ny)\/*
				*/scalar(minority_`y'_5),scalar(me`y'_5_med_ny),scalar(lci`y'_5_med_ny),scalar(uci`y'_5_med_ny)\/*
				*/scalar(minority_`y'_6),scalar(me`y'_6_med_ny),scalar(lci`y'_6_med_ny),scalar(uci`y'_6_med_ny)\/*
				*/scalar(minority_`y'_7),scalar(me`y'_7_med_ny),scalar(lci`y'_7_med_ny),scalar(uci`y'_7_med_ny)\/*
				*/scalar(minority_`y'_8),scalar(me`y'_8_med_ny),scalar(lci`y'_8_med_ny),scalar(uci`y'_8_med_ny)\/*
				*/scalar(minority_`y'_9),scalar(me`y'_9_med_ny),scalar(lci`y'_9_med_ny),scalar(uci`y'_9_med_ny)\/*
				*/scalar(minority_`y'_10),scalar(me`y'_10_med_ny),scalar(lci`y'_10_med_ny),scalar(uci`y'_10_med_ny)\/*
				*/scalar(minority_`y'_11),scalar(me`y'_11_med_ny),scalar(lci`y'_11_med_ny),scalar(uci`y'_11_med_ny))

	matrix colnames me`y'_x_min_med_ny = meminx`y'_med_ny minmex`y'_med_ny minmelcix`y'_med_ny minmeucix`y'_med_ny
	svmat me`y'_x_min_med_ny, names(col)
	
	matrix list me`y'_x_min_med_ny
}

** ordered probit model 

ml model lf myoprobit_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_pminorty_dw_z, nocons) (cut1:) (cut2:) [pw=$wt]
ml maximize
est store min_op_x
	
forvalues y = 0/2 {
	if `y'==0{
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=(-6(1)4)) dydx(f_c9010t_pminorty_dw_z) expression(normal(predict(equation(cut1))-predict(equation(xb1))))
	}
	else if `y'==1{
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=(-6(1)4)) dydx(f_c9010t_pminorty_dw_z) expression((normal(predict(equation(cut2))-predict(equation(xb1))))-(normal(predict(equation(cut1))-predict(equation(xb1)))))
	}
	else if `y'==2{
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=(-6(1)4)) dydx(f_c9010t_pminorty_dw_z) expression(1-(normal(predict(equation(cut2))-predict(equation(xb1)))))
	}
	
	marginsplot, ylabel(-0.2(0.1)0.2, grid) title(" ") xtitle(" ") ytitle(" ") xlabel(-6(2)4) recast(line) recastci(rarea) plotopts(lcolor(black) lwidth(medthick)) ciopts(color(black%50)) ///
	addplot(rarea minmelcix`y'_med_ny minmeucix`y'_med_ny meminx`y'_med_ny, ylabel(-0.2(0.1)0.2, grid) color(black%20)||line minmex`y'_med_ny meminx`y'_med_ny, ylabel(-0.2(0.1)0.2, grid) ///
	lcolor(gs5) lpattern(dash) lwidth(medthick)) legend(off)
	graph export ${pathres}me_y`y'_x_min_med_ny_${date}.pdf, replace
	}

** FIGURE 6 

** W = (hood poverty, hood minority), with covariates

** CT model 

ml model lf myoprobit3_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z, nocons) (cut1:) (lndiff:) /*
						*/(xb2:f_c9010t_perpov_dw_z = x_f_release1 $x_covariates interaction*) (lnvar_v1:) (rho_uv1:) /*
						*/(xb3:f_c9010t_pminorty_dw_z = x_f_release1 $x_covariates interaction*) (lnvar_v2:) (rho_uv2:) (cov_v1v2:) [pw=$wt], /*
						*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
						*/diparm(lnvar_v1, exp label(var_v1)) /*
						*/diparm(lnvar_v2, exp label(var_v2)) 
ml init xb1:x_f_site_bal=0.3039 x_f_site_bos=0.3386 x_f_site_chi=0.1527 x_f_site_la=0.0639 f_c9010t_perpov_dw_z=-0.2998 f_c9010t_pminorty_dw_z=0.3151 cut1:_cons=-0.8165
ml maximize, difficult
est store povmin_t_x

forvalues y = 0/2{
	forvalues i = 1/11{
		scalar poverty_`y'_`i' = `i' - 7
		if `y'==0{
		nlcom -[xb1]f_c9010t_perpov_dw_z * normalden([cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_x*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))
		matrix b = r(b)
		matrix V = r(V)
		scalar povmin1_x_me`y'_`i'_med_ny = b[1,1]
		scalar povmin1_x_v`y'_`i' = sqrt(V[1,1])
		scalar povmin1_x_lci`y'_`i'_med_ny = povmin1_x_me`y'_`i'_med_ny-invnormal(0.975)*povmin1_x_v`y'_`i'
		scalar povmin1_x_uci`y'_`i'_med_ny = povmin1_x_me`y'_`i'_med_ny+invnormal(0.975)*povmin1_x_v`y'_`i'
		}
		else if `y'==1{
		nlcom [xb1]f_c9010t_perpov_dw_z * normalden([cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_x*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex})) - /*
				*/[xb1]f_c9010t_perpov_dw_z * normalden(exp([lndiff]_cons)+[cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_x*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))
		matrix b = r(b)
		matrix V = r(V)
		scalar povmin1_x_me`y'_`i'_med_ny = b[1,1]
		scalar povmin1_x_v`y'_`i' = sqrt(V[1,1])
		scalar povmin1_x_lci`y'_`i'_med_ny = povmin1_x_me`y'_`i'_med_ny-invnormal(0.975)*povmin1_x_v`y'_`i'
		scalar povmin1_x_uci`y'_`i'_med_ny = povmin1_x_me`y'_`i'_med_ny+invnormal(0.975)*povmin1_x_v`y'_`i'
		}
		else if `y'==2{
		nlcom [xb1]f_c9010t_perpov_dw_z * normalden(exp([lndiff]_cons)+[cut1]_cons-poverty_`y'_`i'*[xb1]f_c9010t_perpov_dw_z-f_c9010t_pminorty_dw_z_med_x*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))
		matrix b = r(b)
		matrix V = r(V)
		scalar povmin1_x_me`y'_`i'_med_ny = b[1,1]
		scalar povmin1_x_v`y'_`i' = sqrt(V[1,1])
		scalar povmin1_x_lci`y'_`i'_med_ny = povmin1_x_me`y'_`i'_med_ny-invnormal(0.975)*povmin1_x_v`y'_`i'
		scalar povmin1_x_uci`y'_`i'_med_ny = povmin1_x_me`y'_`i'_med_ny+invnormal(0.975)*povmin1_x_v`y'_`i'
		}
	}

	matrix me`y'_povmin1_x_med_ny = (scalar(poverty_`y'_1),scalar(povmin1_x_me`y'_1_med_ny),scalar(povmin1_x_lci`y'_1_med_ny),scalar(povmin1_x_uci`y'_1_med_ny)\/*
				*/scalar(poverty_`y'_2),scalar(povmin1_x_me`y'_2_med_ny),scalar(povmin1_x_lci`y'_2_med_ny),scalar(povmin1_x_uci`y'_2_med_ny)\/*
				*/scalar(poverty_`y'_3),scalar(povmin1_x_me`y'_3_med_ny),scalar(povmin1_x_lci`y'_3_med_ny),scalar(povmin1_x_uci`y'_3_med_ny)\/*
				*/scalar(poverty_`y'_4),scalar(povmin1_x_me`y'_4_med_ny),scalar(povmin1_x_lci`y'_4_med_ny),scalar(povmin1_x_uci`y'_4_med_ny)\/*
				*/scalar(poverty_`y'_5),scalar(povmin1_x_me`y'_5_med_ny),scalar(povmin1_x_lci`y'_5_med_ny),scalar(povmin1_x_uci`y'_5_med_ny)\/*
				*/scalar(poverty_`y'_6),scalar(povmin1_x_me`y'_6_med_ny),scalar(povmin1_x_lci`y'_6_med_ny),scalar(povmin1_x_uci`y'_6_med_ny)\/*
				*/scalar(poverty_`y'_7),scalar(povmin1_x_me`y'_7_med_ny),scalar(povmin1_x_lci`y'_7_med_ny),scalar(povmin1_x_uci`y'_7_med_ny)\/*
				*/scalar(poverty_`y'_8),scalar(povmin1_x_me`y'_8_med_ny),scalar(povmin1_x_lci`y'_8_med_ny),scalar(povmin1_x_uci`y'_8_med_ny)\/*
				*/scalar(poverty_`y'_9),scalar(povmin1_x_me`y'_9_med_ny),scalar(povmin1_x_lci`y'_9_med_ny),scalar(povmin1_x_uci`y'_9_med_ny)\/*
				*/scalar(poverty_`y'_10),scalar(povmin1_x_me`y'_10_med_ny),scalar(povmin1_x_lci`y'_10_med_ny),scalar(povmin1_x_uci`y'_10_med_ny)\/*
				*/scalar(poverty_`y'_11),scalar(povmin1_x_me`y'_11_med_ny),scalar(povmin1_x_lci`y'_11_med_ny),scalar(povmin1_x_uci`y'_11_med_ny))

	matrix colnames me`y'_povmin1_x_med_ny = povmin1_x_mepov`y'_med_ny povmin1_x_me`y'_med_ny povmin1_x_melci`y'_med_ny povmin1_x_meuci`y'_med_ny
	svmat me`y'_povmin1_x_med_ny, names(col)

	matrix list me`y'_povmin1_x_med_ny
}

forvalues y = 0/2{
	forvalues i = 1/11{
		scalar minority_`y'_`i' = `i' - 7
		if `y'==0{
		nlcom -[xb1]f_c9010t_pminorty_dw_z * normalden([cut1]_cons-f_c9010t_perpov_dw_z_med_x*[xb1]f_c9010t_perpov_dw_z-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))
		matrix b = r(b)
		matrix V = r(V)
		scalar povmin2_x_me`y'_`i'_med_ny = b[1,1]
		scalar povmin2_x_v`y'_`i' = sqrt(V[1,1])
		scalar povmin2_x_lci`y'_`i'_med_ny = povmin2_x_me`y'_`i'_med_ny-invnormal(0.975)*povmin2_x_v`y'_`i'
		scalar povmin2_x_uci`y'_`i'_med_ny = povmin2_x_me`y'_`i'_med_ny+invnormal(0.975)*povmin2_x_v`y'_`i'
		}
		else if `y'==1{
		nlcom [xb1]f_c9010t_pminorty_dw_z * normalden([cut1]_cons-f_c9010t_perpov_dw_z_med_x*[xb1]f_c9010t_perpov_dw_z-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex})) - /*
				*/[xb1]f_c9010t_pminorty_dw_z * normalden(exp([lndiff]_cons)+[cut1]_cons-f_c9010t_perpov_dw_z_med_x*[xb1]f_c9010t_perpov_dw_z-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))
		matrix b = r(b)
		matrix V = r(V)
		scalar povmin2_x_me`y'_`i'_med_ny = b[1,1]
		scalar povmin2_x_v`y'_`i' = sqrt(V[1,1])
		scalar povmin2_x_lci`y'_`i'_med_ny = povmin2_x_me`y'_`i'_med_ny-invnormal(0.975)*povmin2_x_v`y'_`i'
		scalar povmin2_x_uci`y'_`i'_med_ny = povmin2_x_me`y'_`i'_med_ny+invnormal(0.975)*povmin2_x_v`y'_`i'
		}
		else if `y'==2{
		nlcom [xb1]f_c9010t_pminorty_dw_z * normalden(exp([lndiff]_cons)+[cut1]_cons-f_c9010t_perpov_dw_z_med_x*[xb1]f_c9010t_perpov_dw_z-minority_`y'_`i'*[xb1]f_c9010t_pminorty_dw_z-(${xgamma_short}+${xgamma_morex}))
		matrix b = r(b)
		matrix V = r(V)
		scalar povmin2_x_me`y'_`i'_med_ny = b[1,1]
		scalar povmin2_x_v`y'_`i' = sqrt(V[1,1])
		scalar povmin2_x_lci`y'_`i'_med_ny = povmin2_x_me`y'_`i'_med_ny-invnormal(0.975)*povmin2_x_v`y'_`i'
		scalar povmin2_x_uci`y'_`i'_med_ny = povmin2_x_me`y'_`i'_med_ny+invnormal(0.975)*povmin2_x_v`y'_`i'
		}
	}

	matrix me`y'_povmin2_x_med_ny = (scalar(minority_`y'_1),scalar(povmin2_x_me`y'_1_med_ny),scalar(povmin2_x_lci`y'_1_med_ny),scalar(povmin2_x_uci`y'_1_med_ny)\/*
				*/scalar(minority_`y'_2),scalar(povmin2_x_me`y'_2_med_ny),scalar(povmin2_x_lci`y'_2_med_ny),scalar(povmin2_x_uci`y'_2_med_ny)\/*
				*/scalar(minority_`y'_3),scalar(povmin2_x_me`y'_3_med_ny),scalar(povmin2_x_lci`y'_3_med_ny),scalar(povmin2_x_uci`y'_3_med_ny)\/*
				*/scalar(minority_`y'_4),scalar(povmin2_x_me`y'_4_med_ny),scalar(povmin2_x_lci`y'_4_med_ny),scalar(povmin2_x_uci`y'_4_med_ny)\/*
				*/scalar(minority_`y'_5),scalar(povmin2_x_me`y'_5_med_ny),scalar(povmin2_x_lci`y'_5_med_ny),scalar(povmin2_x_uci`y'_5_med_ny)\/*
				*/scalar(minority_`y'_6),scalar(povmin2_x_me`y'_6_med_ny),scalar(povmin2_x_lci`y'_6_med_ny),scalar(povmin2_x_uci`y'_6_med_ny)\/*
				*/scalar(minority_`y'_7),scalar(povmin2_x_me`y'_7_med_ny),scalar(povmin2_x_lci`y'_7_med_ny),scalar(povmin2_x_uci`y'_7_med_ny)\/*
				*/scalar(minority_`y'_8),scalar(povmin2_x_me`y'_8_med_ny),scalar(povmin2_x_lci`y'_8_med_ny),scalar(povmin2_x_uci`y'_8_med_ny)\/*
				*/scalar(minority_`y'_9),scalar(povmin2_x_me`y'_9_med_ny),scalar(povmin2_x_lci`y'_9_med_ny),scalar(povmin2_x_uci`y'_9_med_ny)\/*
				*/scalar(minority_`y'_10),scalar(povmin2_x_me`y'_10_med_ny),scalar(povmin2_x_lci`y'_10_med_ny),scalar(povmin2_x_uci`y'_10_med_ny)\/*
				*/scalar(minority_`y'_11),scalar(povmin2_x_me`y'_11_med_ny),scalar(povmin2_x_lci`y'_11_med_ny),scalar(povmin2_x_uci`y'_11_med_ny))

	matrix colnames me`y'_povmin2_x_med_ny = povmin2_x_memin`y'_med_ny povmin2_x_me`y'_med_ny povmin2_x_melci`y'_med_ny povmin2_x_meuci`y'_med_ny
	svmat me`y'_povmin2_x_med_ny, names(col)

	matrix list me`y'_povmin2_x_med_ny
}

** ordered probit 
	
ml model lf myoprobit_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z, nocons) (cut1:) (cut2:) [pw=$wt]
ml maximize
est store povmin_op_x

forvalues y = 0/2 {
	if `y'==0{
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=.4965963 f_c9010t_perpov_dw_z=(-6(1)4)) dydx(f_c9010t_perpov_dw_z) expression(normal(predict(equation(cut1))-predict(equation(xb1))))
	}
	else if `y'==1{
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=.4965963 f_c9010t_perpov_dw_z=(-6(1)4)) dydx(f_c9010t_perpov_dw_z) expression((normal(predict(equation(cut2))-predict(equation(xb1))))-(normal(predict(equation(cut1))-predict(equation(xb1)))))
	}
	else if `y'==2{
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=.4965963 f_c9010t_perpov_dw_z=(-6(1)4)) dydx(f_c9010t_perpov_dw_z) expression(1-(normal(predict(equation(cut2))-predict(equation(xb1)))))
	}
	
	marginsplot, ylabel(-0.4(0.2)0.4, grid) title(" ") xtitle(" ") ytitle(" ") xlabel(-6(2)4) recast(line) recastci(rarea) plotopts(lcolor(black) lwidth(medthick)) ciopts(color(black%50)) ///
	addplot(rarea povmin1_x_melci`y'_med_ny povmin1_x_meuci`y'_med_ny povmin1_x_mepov`y'_med_ny, ylabel(-0.4(0.2)0.4, grid) color(black%20)||line povmin1_x_me`y'_med_ny ///
	povmin1_x_mepov`y'_med_ny, ylabel(-0.4(0.2)0.4, grid) lcolor(gs5) lpattern(dash) lwidth(medthick)) legend(off)
	graph export ${pathres}me_y`y'_x_povmin1_med_ny_${date}.pdf, replace
	}
	
forvalues y = 0/2 {
	if `y'==0{
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=-.182967 f_c9010t_pminorty_dw_z=(-6(1)4)) dydx(f_c9010t_pminorty_dw_z) expression(normal(predict(equation(cut1))-predict(equation(xb1))))
	}
	else if `y'==1{
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=-.182967 f_c9010t_pminorty_dw_z=(-6(1)4)) dydx(f_c9010t_pminorty_dw_z) expression((normal(predict(equation(cut2))-predict(equation(xb1))))-(normal(predict(equation(cut1))-predict(equation(xb1)))))
	}
	else if `y'==2{
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=-.182967 f_c9010t_pminorty_dw_z=(-6(1)4)) dydx(f_c9010t_pminorty_dw_z) expression(1-(normal(predict(equation(cut2))-predict(equation(xb1)))))
	}
	
	marginsplot, ylabel(-0.4(0.2)0.4, grid) title(" ") xtitle(" ") ytitle(" ") xlabel(-6(2)4) recast(line) recastci(rarea) plotopts(lcolor(black) lwidth(medthick)) ciopts(color(black%50)) ///
	addplot(rarea povmin2_x_melci`y'_med_ny povmin2_x_meuci`y'_med_ny povmin2_x_memin`y'_med_ny, ylabel(-0.4(0.2)0.4, grid) color(black%20)||line povmin2_x_me`y'_med_ny ///
	povmin2_x_memin`y'_med_ny, ylabel(-0.4(0.2)0.4, grid) lcolor(gs5) lpattern(dash) lwidth(medthick)) legend(off)
	graph export ${pathres}me_y`y'_x_povmin2_med_ny_${date}.pdf, replace
	}

scalar drop _all
matrix drop _all
estimates drop _all
