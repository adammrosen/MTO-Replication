
** set seed for bootstrap 

set seed 339487731

** ordered probit parameter estimates, probabilities and marginal effects, CF approach 

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

** W = hood poverty, no covariates

capture noisily program drop my_op_cf_pov_nox 

program my_op_cf_pov_nox, eclass
	tempvar pov_eq2 
	tempvar resid_pov_eq1 resid_pov_eq2 

	regress f_c9010t_perpov_dw_z i.ra_group#i.ra_site x_f_site_bal x_f_site_bos x_f_site_chi x_f_site_la [pw=$wt]
	predict `resid_pov_eq1', residuals 
	
	oprobit happy_scale012 `resid_pov_eq1' f_c9010t_perpov_dw_z x_f_site_bal x_f_site_bos x_f_site_chi x_f_site_la [pw=$wt]
	predict `pov_eq2'
	gen `resid_pov_eq2' = happy_scale012 - `pov_eq2'
	
	correlate `resid_pov_eq1' `resid_pov_eq2'
	local rho_12 = r(rho)
	
	ereturn scalar coef_pov_cf = _b[f_c9010t_perpov_dw_z]*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_bal_cf = _b[x_f_site_bal]*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_bos_cf = _b[x_f_site_bos]*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_chi_cf = _b[x_f_site_chi]*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_la_cf  = _b[x_f_site_la]*((1-((`rho_12')^2))^0.5)
	ereturn scalar cut1		   = _b[/:cut1]*((1-((`rho_12')^2))^0.5)
	ereturn scalar cut2		   = _b[/:cut2]*((1-((`rho_12')^2))^0.5)
	
	local beta_pov  = _b[f_c9010t_perpov_dw_z]*((1-((`rho_12')^2))^0.5)
	local gamma_bal = _b[x_f_site_bal]*((1-((`rho_12')^2))^0.5)
	local gamma_bos = _b[x_f_site_bos]*((1-((`rho_12')^2))^0.5)
	local gamma_chi = _b[x_f_site_chi]*((1-((`rho_12')^2))^0.5)
	local gamma_la  = _b[x_f_site_la]*((1-((`rho_12')^2))^0.5)
	local c1		= _b[/:cut1]*((1-((`rho_12')^2))^0.5)
	local c2		= _b[/:cut2]*((1-((`rho_12')^2))^0.5)
	
end

bootstrap e(coef_pov_cf), reps(100) nodrop: my_op_cf_pov_nox 
bootstrap e(coef_bal_cf), reps(100) nodrop: my_op_cf_pov_nox 
bootstrap e(coef_bos_cf), reps(100) nodrop: my_op_cf_pov_nox 
bootstrap e(coef_chi_cf), reps(100) nodrop: my_op_cf_pov_nox 
bootstrap e(coef_la_cf), reps(100) nodrop: my_op_cf_pov_nox 
bootstrap e(cut1), reps(100) nodrop: my_op_cf_pov_nox 
bootstrap e(cut2), reps(100) nodrop: my_op_cf_pov_nox 

** W = hood poverty, with covariates

capture noisily program drop my_op_cf_pov_x 

program my_op_cf_pov_x, eclass
	tempvar pov_x_eq2 
	tempvar resid_pov_x_eq1 resid_pov_x_eq2 

	regress f_c9010t_perpov_dw_z i.ra_group#i.ra_site $site_covs x_f_release1 $x_covariates [pw=$wt]
	predict `resid_pov_x_eq1', residuals 
	
	oprobit happy_scale012 `resid_pov_x_eq1' f_c9010t_perpov_dw_z $site_covs x_f_release1 $x_covariates [pw=$wt]
	predict `pov_x_eq2'
	gen `resid_pov_x_eq2' = happy_scale012 - `pov_x_eq2'
	
	correlate `resid_pov_x_eq1' `resid_pov_x_eq2'
	local rho_12 = r(rho)
	
	ereturn scalar coef_pov_cf = _b[f_c9010t_perpov_dw_z]*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_bal_cf = _b[x_f_site_bal]	*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_bos_cf = _b[x_f_site_bos]	*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_chi_cf = _b[x_f_site_chi]	*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_la_cf  = _b[x_f_site_la]	*((1-((`rho_12')^2))^0.5)
	ereturn scalar cut1		= _b[/:cut1]*((1-((`rho_12')^2))^0.5)
	ereturn scalar cut2		= _b[/:cut2]*((1-((`rho_12')^2))^0.5)
	
	local beta_pov  = _b[f_c9010t_perpov_dw_z]	*((1-((`rho_12')^2))^0.5)

	local g_bal = _b[x_f_site_bal]	*((1-((`rho_12')^2))^0.5)
	local g_bos = _b[x_f_site_bos]	*((1-((`rho_12')^2))^0.5)
	local g_chi = _b[x_f_site_chi]	*((1-((`rho_12')^2))^0.5)
	local g_la  = _b[x_f_site_la]	*((1-((`rho_12')^2))^0.5)
	
	local g_x_f_release1 				= _b[x_f_release1]				*((1-((`rho_12')^2))^0.5)
    local g_tmp_female 					= _b[tmp_female]				*((1-((`rho_12')^2))^0.5)
    local g_x_rad_ad_le_35 				= _b[x_rad_ad_le_35]			*((1-((`rho_12')^2))^0.5)
    local g_x_rad_ad_36_40 				= _b[x_rad_ad_36_40]			*((1-((`rho_12')^2))^0.5)
    local g_x_rad_ad_41_45 				= _b[x_rad_ad_41_45]			*((1-((`rho_12')^2))^0.5)
    local g_x_rad_ad_46_50 				= _b[x_rad_ad_46_50]			*((1-((`rho_12')^2))^0.5)
	local g_x_rad_ad_ethrace_black_nh 	= _b[x_rad_ad_ethrace_black_nh]	*((1-((`rho_12')^2))^0.5)
    local g_x_rad_ad_ethrace_hisp 		= _b[x_rad_ad_ethrace_hisp]		*((1-((`rho_12')^2))^0.5)
    local g_x_f_ad_nevmarr 				= _b[x_f_ad_nevmarr]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_ad_parentu18 			= _b[x_f_ad_parentu18]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_ad_working 				= _b[x_f_ad_working]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_ad_edinsch 				= _b[x_f_ad_edinsch]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_ad_edgradhs 			= _b[x_f_ad_edgradhs]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_ad_edged 				= _b[x_f_ad_edged]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_afdc 				= _b[x_f_hh_afdc]				*((1-((`rho_12')^2))^0.5)
	local g_rad_svy_bl_totincm_2009d 	= _b[rad_svy_bl_totincm_2009d]	*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_car 					= _b[x_f_hh_car]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_disabl 				= _b[x_f_hh_disabl]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_noteens 				= _b[x_f_hh_noteens]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_size2 				= _b[x_f_hh_size2]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_size3 				= _b[x_f_hh_size3]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_size4 				= _b[x_f_hh_size4]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_victim 				= _b[x_f_hh_victim]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hood_unsafenit 			= _b[x_f_hood_unsafenit]		*((1-((`rho_12')^2))^0.5)
    local g_x_f_hood_verydissat 		= _b[x_f_hood_verydissat]		*((1-((`rho_12')^2))^0.5)
    local g_x_f_hood_5y 				= _b[x_f_hood_5y]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hous_mov3tm 			= _b[x_f_hous_mov3tm]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hood_nofamily 			= _b[x_f_hood_nofamily]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hood_nofriend 			= _b[x_f_hood_nofriend]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hood_chat 				= _b[x_f_hood_chat]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hood_nbrkid 			= _b[x_f_hood_nbrkid]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hous_fndapt 			= _b[x_f_hous_fndapt]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hous_sec8bef 			= _b[x_f_hous_sec8bef]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hous_movdrgs 			= _b[x_f_hous_movdrgs]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hous_movschl 			= _b[x_f_hous_movschl]			*((1-((`rho_12')^2))^0.5)
    local g_cov_hous_movapt 			= _b[cov_hous_movapt]			*((1-((`rho_12')^2))^0.5)
    local g_cov_hous_movjob 			= _b[cov_hous_movjob]			*((1-((`rho_12')^2))^0.5)
	
	local c1		= _b[/:cut1]*((1-((`rho_12')^2))^0.5)
	local c2		= _b[/:cut2]*((1-((`rho_12')^2))^0.5)
	
	local xg_morex 	= (x_f_release1_med_x*`g_x_f_release1') + (tmp_female_med_x*`g_tmp_female') + (x_rad_ad_le_35_med_x*`g_x_rad_ad_le_35') + (x_rad_ad_36_40_med_x*`g_x_rad_ad_36_40') +(x_rad_ad_41_45_med_x*`g_x_rad_ad_41_45') + (x_rad_ad_46_50_med_x*`g_x_rad_ad_46_50') + (x_rad_ad_ethrace_black_nh_med_x*`g_x_rad_ad_ethrace_black_nh') + (x_rad_ad_ethrace_hisp_med_x*`g_x_rad_ad_ethrace_hisp') + (x_f_ad_nevmarr_med_x*`g_x_f_ad_nevmarr') + (x_f_ad_parentu18_med_x*`g_x_f_ad_parentu18') + (x_f_ad_working_med_x*`g_x_f_ad_working') + (x_f_ad_edinsch_med_x*`g_x_f_ad_edinsch') + (x_f_ad_edgradhs_med_x*`g_x_f_ad_edgradhs') + (x_f_ad_edged_med_x*`g_x_f_ad_edged') + (x_f_hh_afdc_med_x*`g_x_f_hh_afdc') + (rad_svy_bl_totincm_2009d_med_x*`g_rad_svy_bl_totincm_2009d') + (x_f_hh_car_med_x*`g_x_f_hh_car') + (x_f_hh_disabl_med_x*`g_x_f_hh_disabl') + (x_f_hh_noteens_med_x*`g_x_f_hh_noteens') + (x_f_hh_size2_med_x*`g_x_f_hh_size2') +(x_f_hh_size3_med_x*`g_x_f_hh_size3') + (x_f_hh_size4_med_x*`g_x_f_hh_size4') + (x_f_hh_victim_med_x*`g_x_f_hh_victim') + (x_f_hood_unsafenit_med_x*`g_x_f_hood_unsafenit') + (x_f_hood_verydissat_med_x*`g_x_f_hood_verydissat') + (x_f_hood_5y_med_x*`g_x_f_hood_5y') + (x_f_hous_mov3tm_med_x*`g_x_f_hous_mov3tm') + (x_f_hood_nofamily_med_x*`g_x_f_hood_nofamily') + (x_f_hood_nofriend_med_x*`g_x_f_hood_nofriend') + (x_f_hood_chat_med_x*`g_x_f_hood_chat') + (x_f_hood_nbrkid_med_x*`g_x_f_hood_nbrkid') + (x_f_hous_fndapt_med_x*`g_x_f_hous_fndapt') + (x_f_hous_sec8bef_med_x*`g_x_f_hous_sec8bef') + (x_f_hous_movdrgs_med_x*`g_x_f_hous_movdrgs') + (x_f_hous_movschl_med_x*`g_x_f_hous_movschl') + (cov_hous_movapt_med_x*`g_cov_hous_movapt') + (cov_hous_movjob_med_x*`g_cov_hous_movjob')
	
	ereturn scalar pr_y0_1  = normal(`c1'-(`beta_pov'*(-6))-`xg_morex')
	ereturn scalar pr_y0_2  = normal(`c1'-(`beta_pov'*(-5))-`xg_morex')
	ereturn scalar pr_y0_3  = normal(`c1'-(`beta_pov'*(-4))-`xg_morex')
	ereturn scalar pr_y0_4  = normal(`c1'-(`beta_pov'*(-3))-`xg_morex')
	ereturn scalar pr_y0_5  = normal(`c1'-(`beta_pov'*(-2))-`xg_morex')
	ereturn scalar pr_y0_6  = normal(`c1'-(`beta_pov'*(-1))-`xg_morex')
	ereturn scalar pr_y0_7  = normal(`c1'-(`beta_pov'*(0))-`xg_morex')
	ereturn scalar pr_y0_8  = normal(`c1'-(`beta_pov'*(1))-`xg_morex')
	ereturn scalar pr_y0_9  = normal(`c1'-(`beta_pov'*(2))-`xg_morex')
	ereturn scalar pr_y0_10 = normal(`c1'-(`beta_pov'*(3))-`xg_morex')
	ereturn scalar pr_y0_11 = normal(`c1'-(`beta_pov'*(4))-`xg_morex')
	
	ereturn scalar pr_y1_1  = normal(`c2'-(`beta_pov'*(-6))-`xg_morex') - normal(`c1'-(`beta_pov'*(-6))-`xg_morex')
	ereturn scalar pr_y1_2  = normal(`c2'-(`beta_pov'*(-5))-`xg_morex') - normal(`c1'-(`beta_pov'*(-5))-`xg_morex')
	ereturn scalar pr_y1_3  = normal(`c2'-(`beta_pov'*(-4))-`xg_morex') - normal(`c1'-(`beta_pov'*(-4))-`xg_morex')
	ereturn scalar pr_y1_4  = normal(`c2'-(`beta_pov'*(-3))-`xg_morex') - normal(`c1'-(`beta_pov'*(-3))-`xg_morex')
	ereturn scalar pr_y1_5  = normal(`c2'-(`beta_pov'*(-2))-`xg_morex') - normal(`c1'-(`beta_pov'*(-2))-`xg_morex')
	ereturn scalar pr_y1_6  = normal(`c2'-(`beta_pov'*(-1))-`xg_morex') - normal(`c1'-(`beta_pov'*(-1))-`xg_morex')
	ereturn scalar pr_y1_7  = normal(`c2'-(`beta_pov'*(0))-`xg_morex') - normal(`c1'-(`beta_pov'*(0))-`xg_morex')
	ereturn scalar pr_y1_8  = normal(`c2'-(`beta_pov'*(1))-`xg_morex') - normal(`c1'-(`beta_pov'*(1))-`xg_morex')
	ereturn scalar pr_y1_9  = normal(`c2'-(`beta_pov'*(2))-`xg_morex') - normal(`c1'-(`beta_pov'*(2))-`xg_morex')
	ereturn scalar pr_y1_10 = normal(`c2'-(`beta_pov'*(3))-`xg_morex') - normal(`c1'-(`beta_pov'*(3))-`xg_morex')
	ereturn scalar pr_y1_11 = normal(`c2'-(`beta_pov'*(4))-`xg_morex') - normal(`c1'-(`beta_pov'*(4))-`xg_morex')
	
	ereturn scalar pr_y2_1  = 1 - normal(`c2'-(`beta_pov'*(-6))-`xg_morex')
	ereturn scalar pr_y2_2  = 1 - normal(`c2'-(`beta_pov'*(-5))-`xg_morex')
	ereturn scalar pr_y2_3  = 1 - normal(`c2'-(`beta_pov'*(-4))-`xg_morex')
	ereturn scalar pr_y2_4  = 1 - normal(`c2'-(`beta_pov'*(-3))-`xg_morex')
	ereturn scalar pr_y2_5  = 1 - normal(`c2'-(`beta_pov'*(-2))-`xg_morex')
	ereturn scalar pr_y2_6  = 1 - normal(`c2'-(`beta_pov'*(-1))-`xg_morex')
	ereturn scalar pr_y2_7  = 1 - normal(`c2'-(`beta_pov'*(0))-`xg_morex')
	ereturn scalar pr_y2_8  = 1 - normal(`c2'-(`beta_pov'*(1))-`xg_morex')
	ereturn scalar pr_y2_9  = 1 - normal(`c2'-(`beta_pov'*(2))-`xg_morex')
	ereturn scalar pr_y2_10 = 1 - normal(`c2'-(`beta_pov'*(3))-`xg_morex')
	ereturn scalar pr_y2_11 = 1 - normal(`c2'-(`beta_pov'*(4))-`xg_morex')
	
	ereturn scalar me_y0_1  = -(`beta_pov')*(normalden(`c1'-(`beta_pov'*(-6))-`xg_morex'))
	ereturn scalar me_y0_2  = -(`beta_pov')*(normalden(`c1'-(`beta_pov'*(-5))-`xg_morex'))
	ereturn scalar me_y0_3  = -(`beta_pov')*(normalden(`c1'-(`beta_pov'*(-4))-`xg_morex'))
	ereturn scalar me_y0_4  = -(`beta_pov')*(normalden(`c1'-(`beta_pov'*(-3))-`xg_morex'))
	ereturn scalar me_y0_5  = -(`beta_pov')*(normalden(`c1'-(`beta_pov'*(-2))-`xg_morex'))
	ereturn scalar me_y0_6  = -(`beta_pov')*(normalden(`c1'-(`beta_pov'*(-1))-`xg_morex'))
	ereturn scalar me_y0_7  = -(`beta_pov')*(normalden(`c1'-(`beta_pov'*(0))-`xg_morex'))
	ereturn scalar me_y0_8  = -(`beta_pov')*(normalden(`c1'-(`beta_pov'*(1))-`xg_morex'))
	ereturn scalar me_y0_9  = -(`beta_pov')*(normalden(`c1'-(`beta_pov'*(2))-`xg_morex'))
	ereturn scalar me_y0_10 = -(`beta_pov')*(normalden(`c1'-(`beta_pov'*(3))-`xg_morex'))
	ereturn scalar me_y0_11 = -(`beta_pov')*(normalden(`c1'-(`beta_pov'*(4))-`xg_morex'))
	
	ereturn scalar me_y1_1  = -(`beta_pov')*(normalden(`c2'-(`beta_pov'*(-6))-`xg_morex') - normalden(`c1'-(`beta_pov'*(-6))-`xg_morex'))
	ereturn scalar me_y1_2  = -(`beta_pov')*(normalden(`c2'-(`beta_pov'*(-5))-`xg_morex') - normalden(`c1'-(`beta_pov'*(-5))-`xg_morex'))
	ereturn scalar me_y1_3  = -(`beta_pov')*(normalden(`c2'-(`beta_pov'*(-4))-`xg_morex') - normalden(`c1'-(`beta_pov'*(-4))-`xg_morex'))
	ereturn scalar me_y1_4  = -(`beta_pov')*(normalden(`c2'-(`beta_pov'*(-3))-`xg_morex') - normalden(`c1'-(`beta_pov'*(-3))-`xg_morex'))
	ereturn scalar me_y1_5  = -(`beta_pov')*(normalden(`c2'-(`beta_pov'*(-2))-`xg_morex') - normalden(`c1'-(`beta_pov'*(-2))-`xg_morex'))
	ereturn scalar me_y1_6  = -(`beta_pov')*(normalden(`c2'-(`beta_pov'*(-1))-`xg_morex') - normalden(`c1'-(`beta_pov'*(-1))-`xg_morex'))
	ereturn scalar me_y1_7  = -(`beta_pov')*(normalden(`c2'-(`beta_pov'*(0))-`xg_morex') - normalden(`c1'-(`beta_pov'*(0))-`xg_morex'))
	ereturn scalar me_y1_8  = -(`beta_pov')*(normalden(`c2'-(`beta_pov'*(1))-`xg_morex') - normalden(`c1'-(`beta_pov'*(1))-`xg_morex'))
	ereturn scalar me_y1_9  = -(`beta_pov')*(normalden(`c2'-(`beta_pov'*(2))-`xg_morex') - normalden(`c1'-(`beta_pov'*(2))-`xg_morex'))
	ereturn scalar me_y1_10 = -(`beta_pov')*(normalden(`c2'-(`beta_pov'*(3))-`xg_morex') - normalden(`c1'-(`beta_pov'*(3))-`xg_morex'))
	ereturn scalar me_y1_11 = -(`beta_pov')*(normalden(`c2'-(`beta_pov'*(4))-`xg_morex') - normalden(`c1'-(`beta_pov'*(4))-`xg_morex'))
	
	ereturn scalar me_y2_1  = -(`beta_pov')*(-normalden(`c2'-(`beta_pov'*(-6))-`xg_morex'))
	ereturn scalar me_y2_2  = -(`beta_pov')*(-normalden(`c2'-(`beta_pov'*(-5))-`xg_morex'))
	ereturn scalar me_y2_3  = -(`beta_pov')*(-normalden(`c2'-(`beta_pov'*(-4))-`xg_morex'))
	ereturn scalar me_y2_4  = -(`beta_pov')*(-normalden(`c2'-(`beta_pov'*(-3))-`xg_morex'))
	ereturn scalar me_y2_5  = -(`beta_pov')*(-normalden(`c2'-(`beta_pov'*(-2))-`xg_morex'))
	ereturn scalar me_y2_6  = -(`beta_pov')*(-normalden(`c2'-(`beta_pov'*(-1))-`xg_morex'))
	ereturn scalar me_y2_7  = -(`beta_pov')*(-normalden(`c2'-(`beta_pov'*(0))-`xg_morex'))
	ereturn scalar me_y2_8  = -(`beta_pov')*(-normalden(`c2'-(`beta_pov'*(1))-`xg_morex'))
	ereturn scalar me_y2_9  = -(`beta_pov')*(-normalden(`c2'-(`beta_pov'*(2))-`xg_morex'))
	ereturn scalar me_y2_10 = -(`beta_pov')*(-normalden(`c2'-(`beta_pov'*(3))-`xg_morex'))
	ereturn scalar me_y2_11 = -(`beta_pov')*(-normalden(`c2'-(`beta_pov'*(4))-`xg_morex'))
	
end

bootstrap e(coef_pov_cf), reps(100) nodrop: my_op_cf_pov_x 
bootstrap e(coef_bal_cf), reps(100) nodrop: my_op_cf_pov_x 
bootstrap e(coef_bos_cf), reps(100) nodrop: my_op_cf_pov_x 
bootstrap e(coef_chi_cf), reps(100) nodrop: my_op_cf_pov_x 
bootstrap e(coef_la_cf), reps(100) nodrop: my_op_cf_pov_x 
bootstrap e(cut1), reps(100) nodrop: my_op_cf_pov_x 
bootstrap e(cut2), reps(100) nodrop: my_op_cf_pov_x 

*** conditional probabilities that Y = 0,1,2

forvalues y = 0(1)2{ 
	forvalues i = 1(1)11{
    
		scalar pov_pr_y`y'_`i' = `i'-7

		bootstrap e(pr_y`y'_`i'), reps(100) nodrop: my_op_cf_pov_x  
		scalar pr_y`y'_`i' 		= e(b)[1,1]
		scalar pr_y`y'_`i'_v	= e(V)[1,1]
	
		scalar pr_y`y'_`i'_lci = pr_y`y'_`i'-invnormal(0.975)*(sqrt(pr_y`y'_`i'_v))
		scalar pr_y`y'_`i'_uci = pr_y`y'_`i'+invnormal(0.975)*(sqrt(pr_y`y'_`i'_v))

	}

	matrix pov_pr_y`y'_x = (scalar(pov_pr_y`y'_1),scalar(pr_y`y'_1),scalar(pr_y`y'_1_lci),scalar(pr_y`y'_1_uci)\/*
				*/scalar(pov_pr_y`y'_2),scalar(pr_y`y'_2),scalar(pr_y`y'_2_lci),scalar(pr_y`y'_2_uci)\/*
				*/scalar(pov_pr_y`y'_3),scalar(pr_y`y'_3),scalar(pr_y`y'_3_lci),scalar(pr_y`y'_3_uci)\/*
				*/scalar(pov_pr_y`y'_4),scalar(pr_y`y'_4),scalar(pr_y`y'_4_lci),scalar(pr_y`y'_4_uci)\/*
				*/scalar(pov_pr_y`y'_5),scalar(pr_y`y'_5),scalar(pr_y`y'_5_lci),scalar(pr_y`y'_5_uci)\/*
				*/scalar(pov_pr_y`y'_6),scalar(pr_y`y'_6),scalar(pr_y`y'_6_lci),scalar(pr_y`y'_6_uci)\/*
				*/scalar(pov_pr_y`y'_7),scalar(pr_y`y'_7),scalar(pr_y`y'_7_lci),scalar(pr_y`y'_7_uci)\/*
				*/scalar(pov_pr_y`y'_8),scalar(pr_y`y'_8),scalar(pr_y`y'_8_lci),scalar(pr_y`y'_8_uci)\/*
				*/scalar(pov_pr_y`y'_9),scalar(pr_y`y'_9),scalar(pr_y`y'_9_lci),scalar(pr_y`y'_9_uci)\/*
				*/scalar(pov_pr_y`y'_10),scalar(pr_y`y'_10),scalar(pr_y`y'_10_lci),scalar(pr_y`y'_10_uci)\/*
				*/scalar(pov_pr_y`y'_11),scalar(pr_y`y'_11),scalar(pr_y`y'_11_lci),scalar(pr_y`y'_11_uci))

	matrix colnames pov_pr_y`y'_x = poverty_pr_y`y'_x pov_probs_y`y'_x pov_pr_lc_y`y'_x pov_pr_uc_y`y'_x
	svmat pov_pr_y`y'_x, names(col)
	twoway rcap pov_pr_lc_y`y'_x pov_pr_uc_y`y'_x poverty_pr_y`y'_x, lcolor(black)||scatter pov_probs_y`y'_x poverty_pr_y`y'_x, mcolor(black)||line pov_probs_y`y'_x poverty_pr_y`y'_x, /*
	*/lcolor(black) title(" ") xtitle(" ") ytitle(" ") ylabel(-0.4(0.4)1.2) legend(off)
	
	graph export ${pathres}prob_CF_op_y`y'_x_pov_med_ny_${date}.eps, replace
	
	scalar drop pov_pr_y`y'_1 pov_pr_y`y'_2 pov_pr_y`y'_3 pov_pr_y`y'_4 pov_pr_y`y'_5 pov_pr_y`y'_6 pov_pr_y`y'_7 pov_pr_y`y'_8 pov_pr_y`y'_9 pov_pr_y`y'_10 pov_pr_y`y'_11 /*
			*/ pr_y`y'_1 pr_y`y'_2 pr_y`y'_3 pr_y`y'_4 pr_y`y'_5 pr_y`y'_6 pr_y`y'_7 pr_y`y'_8 pr_y`y'_9 pr_y`y'_10 pr_y`y'_11 /*
			*/ pr_y`y'_1_lci pr_y`y'_2_lci pr_y`y'_3_lci pr_y`y'_4_lci pr_y`y'_5_lci pr_y`y'_6_lci pr_y`y'_7_lci pr_y`y'_8_lci pr_y`y'_9_lci pr_y`y'_10_lci pr_y`y'_11_lci /*
			*/ pr_y`y'_1_uci pr_y`y'_2_uci pr_y`y'_3_uci pr_y`y'_4_uci pr_y`y'_5_uci pr_y`y'_6_uci pr_y`y'_7_uci pr_y`y'_8_uci pr_y`y'_9_uci pr_y`y'_10_uci pr_y`y'_11_uci 
}

*** conditional marginal effects for Y = 0,1,2

forvalues y = 0(1)2{ 
	forvalues i = 1(1)11{
    
		scalar pov_me_y`y'_`i' = `i'-7

		bootstrap e(me_y`y'_`i'), reps(100) nodrop: my_op_cf_pov_x 
		scalar me_y`y'_`i' 		= e(b)[1,1]
		scalar me_y`y'_`i'_v	= e(V)[1,1]
	
		scalar me_y`y'_`i'_lci = me_y`y'_`i'-invnormal(0.975)*(sqrt(me_y`y'_`i'_v))
		scalar me_y`y'_`i'_uci = me_y`y'_`i'+invnormal(0.975)*(sqrt(me_y`y'_`i'_v))

	}

	matrix pov_me_y`y'_x = (scalar(pov_me_y`y'_1),scalar(me_y`y'_1),scalar(me_y`y'_1_lci),scalar(me_y`y'_1_uci)\/*
				*/scalar(pov_me_y`y'_2),scalar(me_y`y'_2),scalar(me_y`y'_2_lci),scalar(me_y`y'_2_uci)\/*
				*/scalar(pov_me_y`y'_3),scalar(me_y`y'_3),scalar(me_y`y'_3_lci),scalar(me_y`y'_3_uci)\/*
				*/scalar(pov_me_y`y'_4),scalar(me_y`y'_4),scalar(me_y`y'_4_lci),scalar(me_y`y'_4_uci)\/*
				*/scalar(pov_me_y`y'_5),scalar(me_y`y'_5),scalar(me_y`y'_5_lci),scalar(me_y`y'_5_uci)\/*
				*/scalar(pov_me_y`y'_6),scalar(me_y`y'_6),scalar(me_y`y'_6_lci),scalar(me_y`y'_6_uci)\/*
				*/scalar(pov_me_y`y'_7),scalar(me_y`y'_7),scalar(me_y`y'_7_lci),scalar(me_y`y'_7_uci)\/*
				*/scalar(pov_me_y`y'_8),scalar(me_y`y'_8),scalar(me_y`y'_8_lci),scalar(me_y`y'_8_uci)\/*
				*/scalar(pov_me_y`y'_9),scalar(me_y`y'_9),scalar(me_y`y'_9_lci),scalar(me_y`y'_9_uci)\/*
				*/scalar(pov_me_y`y'_10),scalar(me_y`y'_10),scalar(me_y`y'_10_lci),scalar(me_y`y'_10_uci)\/*
				*/scalar(pov_me_y`y'_11),scalar(me_y`y'_11),scalar(me_y`y'_11_lci),scalar(me_y`y'_11_uci))

	matrix colnames pov_me_y`y'_x = poverty_me_y`y'_x pov_meffects_y`y'_x pov_me_lc_y`y'_x pov_me_uc_y`y'_x
	svmat pov_me_y`y'_x, names(col)
	twoway rcap pov_me_lc_y`y'_x pov_me_uc_y`y'_x poverty_me_y`y'_x, lcolor(black)||scatter pov_meffects_y`y'_x poverty_me_y`y'_x, mcolor(black)||line pov_meffects_y`y'_x poverty_me_y`y'_x, /*
	*/lcolor(black) title(" ") xtitle(" ") ytitle(" ") ylabel(-0.2(0.1)0.2) legend(off)
	
	graph export ${pathres}me_CF_op_y`y'_x_pov_med_ny_${date}.eps, replace
	
	scalar drop pov_me_y`y'_1 pov_me_y`y'_2 pov_me_y`y'_3 pov_me_y`y'_4 pov_me_y`y'_5 pov_me_y`y'_6 pov_me_y`y'_7 pov_me_y`y'_8 pov_me_y`y'_9 pov_me_y`y'_10 pov_me_y`y'_11 /*
			*/ me_y`y'_1 me_y`y'_2 me_y`y'_3 me_y`y'_4 me_y`y'_5 me_y`y'_6 me_y`y'_7 me_y`y'_8 me_y`y'_9 me_y`y'_10 me_y`y'_11 /*
			*/ me_y`y'_1_lci me_y`y'_2_lci me_y`y'_3_lci me_y`y'_4_lci me_y`y'_5_lci me_y`y'_6_lci me_y`y'_7_lci me_y`y'_8_lci me_y`y'_9_lci me_y`y'_10_lci me_y`y'_11_lci /*
			*/ me_y`y'_1_uci me_y`y'_2_uci me_y`y'_3_uci me_y`y'_4_uci me_y`y'_5_uci me_y`y'_6_uci me_y`y'_7_uci me_y`y'_8_uci me_y`y'_9_uci me_y`y'_10_uci me_y`y'_11_uci 
}

** W = hood minority, no covariates

capture noisily program drop my_op_cf_min_nox 

program my_op_cf_min_nox, eclass
	tempvar min_eq2 
	tempvar resid_min_eq1 resid_min_eq2 

	regress f_c9010t_pminorty_dw_z i.ra_group#i.ra_site x_f_site_bal x_f_site_bos x_f_site_chi x_f_site_la [pw=$wt]
	predict `resid_min_eq1', residuals 
	
	oprobit happy_scale012 `resid_min_eq1' f_c9010t_pminorty_dw_z x_f_site_bal x_f_site_bos x_f_site_chi x_f_site_la [pw=$wt]
	predict `min_eq2'
	gen `resid_min_eq2' = happy_scale012 - `min_eq2'
	
	correlate `resid_min_eq1' `resid_min_eq2'
	local rho_12 = r(rho)
	
	ereturn scalar coef_min_cf = _b[f_c9010t_pminorty_dw_z]*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_bal_cf = _b[x_f_site_bal]*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_bos_cf = _b[x_f_site_bos]*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_chi_cf = _b[x_f_site_chi]*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_la_cf  = _b[x_f_site_la]*((1-((`rho_12')^2))^0.5)
	ereturn scalar cut1		= _b[/:cut1]*((1-((`rho_12')^2))^0.5)
	ereturn scalar cut2		= _b[/:cut2]*((1-((`rho_12')^2))^0.5)
	
	local beta_min  = _b[f_c9010t_pminorty_dw_z]*((1-((`rho_12')^2))^0.5)
	local gamma_bal = _b[x_f_site_bal]*((1-((`rho_12')^2))^0.5)
	local gamma_bos = _b[x_f_site_bos]*((1-((`rho_12')^2))^0.5)
	local gamma_chi = _b[x_f_site_chi]*((1-((`rho_12')^2))^0.5)
	local gamma_la  = _b[x_f_site_la]*((1-((`rho_12')^2))^0.5)
	
	local c1		= _b[/:cut1]*((1-((`rho_12')^2))^0.5)
	local c2		= _b[/:cut2]*((1-((`rho_12')^2))^0.5)
	
end

bootstrap e(coef_min_cf), reps(100) nodrop: my_op_cf_min_nox 
bootstrap e(coef_bal_cf), reps(100) nodrop: my_op_cf_min_nox 
bootstrap e(coef_bos_cf), reps(100) nodrop: my_op_cf_min_nox 
bootstrap e(coef_chi_cf), reps(100) nodrop: my_op_cf_min_nox 
bootstrap e(coef_la_cf), reps(100) nodrop: my_op_cf_min_nox 
bootstrap e(cut1), reps(100) nodrop: my_op_cf_min_nox 
bootstrap e(cut2), reps(100) nodrop: my_op_cf_min_nox 

** W = hood minority, with covariates

capture noisily program drop my_op_cf_min_x 

program my_op_cf_min_x, eclass
	tempvar min_x_eq2 
	tempvar resid_min_x_eq1 resid_min_x_eq2 

	regress f_c9010t_pminorty_dw_z i.ra_group#i.ra_site $site_covs x_f_release1 $x_covariates [pw=$wt]
	predict `resid_min_x_eq1', residuals 
	
	oprobit happy_scale012 `resid_min_x_eq1' f_c9010t_pminorty_dw_z $site_covs x_f_release1 $x_covariates [pw=$wt]
	predict `min_x_eq2'
	gen `resid_min_x_eq2' = happy_scale012 - `min_x_eq2'
	
	correlate `resid_min_x_eq1' `resid_min_x_eq2'
	local rho_12 = r(rho)
	
	ereturn scalar coef_min_cf = _b[f_c9010t_pminorty_dw_z]*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_bal_cf = _b[x_f_site_bal]	*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_bos_cf = _b[x_f_site_bos]	*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_chi_cf = _b[x_f_site_chi]	*((1-((`rho_12')^2))^0.5)
	ereturn scalar coef_la_cf  = _b[x_f_site_la]	*((1-((`rho_12')^2))^0.5)
	ereturn scalar cut1		= _b[/:cut1]*((1-((`rho_12')^2))^0.5)
	ereturn scalar cut2		= _b[/:cut2]*((1-((`rho_12')^2))^0.5)
	
	local beta_min  = _b[f_c9010t_pminorty_dw_z]	*((1-((`rho_12')^2))^0.5)

	local g_bal = _b[x_f_site_bal]	*((1-((`rho_12')^2))^0.5)
	local g_bos = _b[x_f_site_bos]	*((1-((`rho_12')^2))^0.5)
	local g_chi = _b[x_f_site_chi]	*((1-((`rho_12')^2))^0.5)
	local g_la  = _b[x_f_site_la]	*((1-((`rho_12')^2))^0.5)
	
	local g_x_f_release1 				= _b[x_f_release1]				*((1-((`rho_12')^2))^0.5)
    local g_tmp_female 					= _b[tmp_female]				*((1-((`rho_12')^2))^0.5)
    local g_x_rad_ad_le_35 				= _b[x_rad_ad_le_35]			*((1-((`rho_12')^2))^0.5)
    local g_x_rad_ad_36_40 				= _b[x_rad_ad_36_40]			*((1-((`rho_12')^2))^0.5)
    local g_x_rad_ad_41_45 				= _b[x_rad_ad_41_45]			*((1-((`rho_12')^2))^0.5)
    local g_x_rad_ad_46_50 				= _b[x_rad_ad_46_50]			*((1-((`rho_12')^2))^0.5)
	local g_x_rad_ad_ethrace_black_nh 	= _b[x_rad_ad_ethrace_black_nh]	*((1-((`rho_12')^2))^0.5)
    local g_x_rad_ad_ethrace_hisp 		= _b[x_rad_ad_ethrace_hisp]		*((1-((`rho_12')^2))^0.5)
    local g_x_f_ad_nevmarr 				= _b[x_f_ad_nevmarr]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_ad_parentu18 			= _b[x_f_ad_parentu18]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_ad_working 				= _b[x_f_ad_working]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_ad_edinsch 				= _b[x_f_ad_edinsch]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_ad_edgradhs 			= _b[x_f_ad_edgradhs]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_ad_edged 				= _b[x_f_ad_edged]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_afdc 				= _b[x_f_hh_afdc]				*((1-((`rho_12')^2))^0.5)
	local g_rad_svy_bl_totincm_2009d 	= _b[rad_svy_bl_totincm_2009d]	*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_car 					= _b[x_f_hh_car]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_disabl 				= _b[x_f_hh_disabl]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_noteens 				= _b[x_f_hh_noteens]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_size2 				= _b[x_f_hh_size2]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_size3 				= _b[x_f_hh_size3]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_size4 				= _b[x_f_hh_size4]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hh_victim 				= _b[x_f_hh_victim]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hood_unsafenit 			= _b[x_f_hood_unsafenit]		*((1-((`rho_12')^2))^0.5)
    local g_x_f_hood_verydissat 		= _b[x_f_hood_verydissat]		*((1-((`rho_12')^2))^0.5)
    local g_x_f_hood_5y 				= _b[x_f_hood_5y]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hous_mov3tm 			= _b[x_f_hous_mov3tm]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hood_nofamily 			= _b[x_f_hood_nofamily]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hood_nofriend 			= _b[x_f_hood_nofriend]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hood_chat 				= _b[x_f_hood_chat]				*((1-((`rho_12')^2))^0.5)
    local g_x_f_hood_nbrkid 			= _b[x_f_hood_nbrkid]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hous_fndapt 			= _b[x_f_hous_fndapt]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hous_sec8bef 			= _b[x_f_hous_sec8bef]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hous_movdrgs 			= _b[x_f_hous_movdrgs]			*((1-((`rho_12')^2))^0.5)
    local g_x_f_hous_movschl 			= _b[x_f_hous_movschl]			*((1-((`rho_12')^2))^0.5)
    local g_cov_hous_movapt 			= _b[cov_hous_movapt]			*((1-((`rho_12')^2))^0.5)
    local g_cov_hous_movjob 			= _b[cov_hous_movjob]			*((1-((`rho_12')^2))^0.5)
	
	local c1		= _b[/:cut1]*((1-((`rho_12')^2))^0.5)
	local c2		= _b[/:cut2]*((1-((`rho_12')^2))^0.5)
	
	local xg_morex 	= (x_f_release1_med_x*`g_x_f_release1') + (tmp_female_med_x*`g_tmp_female') + (x_rad_ad_le_35_med_x*`g_x_rad_ad_le_35') + (x_rad_ad_36_40_med_x*`g_x_rad_ad_36_40') +(x_rad_ad_41_45_med_x*`g_x_rad_ad_41_45') + (x_rad_ad_46_50_med_x*`g_x_rad_ad_46_50') + (x_rad_ad_ethrace_black_nh_med_x*`g_x_rad_ad_ethrace_black_nh') + (x_rad_ad_ethrace_hisp_med_x*`g_x_rad_ad_ethrace_hisp') + (x_f_ad_nevmarr_med_x*`g_x_f_ad_nevmarr') + (x_f_ad_parentu18_med_x*`g_x_f_ad_parentu18') + (x_f_ad_working_med_x*`g_x_f_ad_working') + (x_f_ad_edinsch_med_x*`g_x_f_ad_edinsch') + (x_f_ad_edgradhs_med_x*`g_x_f_ad_edgradhs') + (x_f_ad_edged_med_x*`g_x_f_ad_edged') + (x_f_hh_afdc_med_x*`g_x_f_hh_afdc') + (rad_svy_bl_totincm_2009d_med_x*`g_rad_svy_bl_totincm_2009d') + (x_f_hh_car_med_x*`g_x_f_hh_car') + (x_f_hh_disabl_med_x*`g_x_f_hh_disabl') + (x_f_hh_noteens_med_x*`g_x_f_hh_noteens') + (x_f_hh_size2_med_x*`g_x_f_hh_size2') +(x_f_hh_size3_med_x*`g_x_f_hh_size3') + (x_f_hh_size4_med_x*`g_x_f_hh_size4') + (x_f_hh_victim_med_x*`g_x_f_hh_victim') + (x_f_hood_unsafenit_med_x*`g_x_f_hood_unsafenit') + (x_f_hood_verydissat_med_x*`g_x_f_hood_verydissat') + (x_f_hood_5y_med_x*`g_x_f_hood_5y') + (x_f_hous_mov3tm_med_x*`g_x_f_hous_mov3tm') + (x_f_hood_nofamily_med_x*`g_x_f_hood_nofamily') + (x_f_hood_nofriend_med_x*`g_x_f_hood_nofriend') + (x_f_hood_chat_med_x*`g_x_f_hood_chat') + (x_f_hood_nbrkid_med_x*`g_x_f_hood_nbrkid') + (x_f_hous_fndapt_med_x*`g_x_f_hous_fndapt') + (x_f_hous_sec8bef_med_x*`g_x_f_hous_sec8bef') + (x_f_hous_movdrgs_med_x*`g_x_f_hous_movdrgs') + (x_f_hous_movschl_med_x*`g_x_f_hous_movschl') + (cov_hous_movapt_med_x*`g_cov_hous_movapt') + (cov_hous_movjob_med_x*`g_cov_hous_movjob')
	
	ereturn scalar pr_y0_1  = normal(`c1'-(`beta_min'*(-6))-`xg_morex')
	ereturn scalar pr_y0_2  = normal(`c1'-(`beta_min'*(-5))-`xg_morex')
	ereturn scalar pr_y0_3  = normal(`c1'-(`beta_min'*(-4))-`xg_morex')
	ereturn scalar pr_y0_4  = normal(`c1'-(`beta_min'*(-3))-`xg_morex')
	ereturn scalar pr_y0_5  = normal(`c1'-(`beta_min'*(-2))-`xg_morex')
	ereturn scalar pr_y0_6  = normal(`c1'-(`beta_min'*(-1))-`xg_morex')
	ereturn scalar pr_y0_7  = normal(`c1'-(`beta_min'*(0))-`xg_morex')
	ereturn scalar pr_y0_8  = normal(`c1'-(`beta_min'*(1))-`xg_morex')
	ereturn scalar pr_y0_9  = normal(`c1'-(`beta_min'*(2))-`xg_morex')
	ereturn scalar pr_y0_10 = normal(`c1'-(`beta_min'*(3))-`xg_morex')
	ereturn scalar pr_y0_11 = normal(`c1'-(`beta_min'*(4))-`xg_morex')
	
	ereturn scalar pr_y1_1  = normal(`c2'-(`beta_min'*(-6))-`xg_morex') - normal(`c1'-(`beta_min'*(-6))-`xg_morex')
	ereturn scalar pr_y1_2  = normal(`c2'-(`beta_min'*(-5))-`xg_morex') - normal(`c1'-(`beta_min'*(-5))-`xg_morex')
	ereturn scalar pr_y1_3  = normal(`c2'-(`beta_min'*(-4))-`xg_morex') - normal(`c1'-(`beta_min'*(-4))-`xg_morex')
	ereturn scalar pr_y1_4  = normal(`c2'-(`beta_min'*(-3))-`xg_morex') - normal(`c1'-(`beta_min'*(-3))-`xg_morex')
	ereturn scalar pr_y1_5  = normal(`c2'-(`beta_min'*(-2))-`xg_morex') - normal(`c1'-(`beta_min'*(-2))-`xg_morex')
	ereturn scalar pr_y1_6  = normal(`c2'-(`beta_min'*(-1))-`xg_morex') - normal(`c1'-(`beta_min'*(-1))-`xg_morex')
	ereturn scalar pr_y1_7  = normal(`c2'-(`beta_min'*(0))-`xg_morex') - normal(`c1'-(`beta_min'*(0))-`xg_morex')
	ereturn scalar pr_y1_8  = normal(`c2'-(`beta_min'*(1))-`xg_morex') - normal(`c1'-(`beta_min'*(1))-`xg_morex')
	ereturn scalar pr_y1_9  = normal(`c2'-(`beta_min'*(2))-`xg_morex') - normal(`c1'-(`beta_min'*(2))-`xg_morex')
	ereturn scalar pr_y1_10 = normal(`c2'-(`beta_min'*(3))-`xg_morex') - normal(`c1'-(`beta_min'*(3))-`xg_morex')
	ereturn scalar pr_y1_11 = normal(`c2'-(`beta_min'*(4))-`xg_morex') - normal(`c1'-(`beta_min'*(4))-`xg_morex')
	
	ereturn scalar pr_y2_1  = 1 - normal(`c2'-(`beta_min'*(-6))-`xg_morex')
	ereturn scalar pr_y2_2  = 1 - normal(`c2'-(`beta_min'*(-5))-`xg_morex')
	ereturn scalar pr_y2_3  = 1 - normal(`c2'-(`beta_min'*(-4))-`xg_morex')
	ereturn scalar pr_y2_4  = 1 - normal(`c2'-(`beta_min'*(-3))-`xg_morex')
	ereturn scalar pr_y2_5  = 1 - normal(`c2'-(`beta_min'*(-2))-`xg_morex')
	ereturn scalar pr_y2_6  = 1 - normal(`c2'-(`beta_min'*(-1))-`xg_morex')
	ereturn scalar pr_y2_7  = 1 - normal(`c2'-(`beta_min'*(0))-`xg_morex')
	ereturn scalar pr_y2_8  = 1 - normal(`c2'-(`beta_min'*(1))-`xg_morex')
	ereturn scalar pr_y2_9  = 1 - normal(`c2'-(`beta_min'*(2))-`xg_morex')
	ereturn scalar pr_y2_10 = 1 - normal(`c2'-(`beta_min'*(3))-`xg_morex')
	ereturn scalar pr_y2_11 = 1 - normal(`c2'-(`beta_min'*(4))-`xg_morex')
	
	ereturn scalar me_y0_1  = -(`beta_min')*(normalden(`c1'-(`beta_min'*(-6))-`xg_morex'))
	ereturn scalar me_y0_2  = -(`beta_min')*(normalden(`c1'-(`beta_min'*(-5))-`xg_morex'))
	ereturn scalar me_y0_3  = -(`beta_min')*(normalden(`c1'-(`beta_min'*(-4))-`xg_morex'))
	ereturn scalar me_y0_4  = -(`beta_min')*(normalden(`c1'-(`beta_min'*(-3))-`xg_morex'))
	ereturn scalar me_y0_5  = -(`beta_min')*(normalden(`c1'-(`beta_min'*(-2))-`xg_morex'))
	ereturn scalar me_y0_6  = -(`beta_min')*(normalden(`c1'-(`beta_min'*(-1))-`xg_morex'))
	ereturn scalar me_y0_7  = -(`beta_min')*(normalden(`c1'-(`beta_min'*(0))-`xg_morex'))
	ereturn scalar me_y0_8  = -(`beta_min')*(normalden(`c1'-(`beta_min'*(1))-`xg_morex'))
	ereturn scalar me_y0_9  = -(`beta_min')*(normalden(`c1'-(`beta_min'*(2))-`xg_morex'))
	ereturn scalar me_y0_10 = -(`beta_min')*(normalden(`c1'-(`beta_min'*(3))-`xg_morex'))
	ereturn scalar me_y0_11 = -(`beta_min')*(normalden(`c1'-(`beta_min'*(4))-`xg_morex'))
	
	ereturn scalar me_y1_1  = -(`beta_min')*(normalden(`c2'-(`beta_min'*(-6))-`xg_morex') - normalden(`c1'-(`beta_min'*(-6))-`xg_morex'))
	ereturn scalar me_y1_2  = -(`beta_min')*(normalden(`c2'-(`beta_min'*(-5))-`xg_morex') - normalden(`c1'-(`beta_min'*(-5))-`xg_morex'))
	ereturn scalar me_y1_3  = -(`beta_min')*(normalden(`c2'-(`beta_min'*(-4))-`xg_morex') - normalden(`c1'-(`beta_min'*(-4))-`xg_morex'))
	ereturn scalar me_y1_4  = -(`beta_min')*(normalden(`c2'-(`beta_min'*(-3))-`xg_morex') - normalden(`c1'-(`beta_min'*(-3))-`xg_morex'))
	ereturn scalar me_y1_5  = -(`beta_min')*(normalden(`c2'-(`beta_min'*(-2))-`xg_morex') - normalden(`c1'-(`beta_min'*(-2))-`xg_morex'))
	ereturn scalar me_y1_6  = -(`beta_min')*(normalden(`c2'-(`beta_min'*(-1))-`xg_morex') - normalden(`c1'-(`beta_min'*(-1))-`xg_morex'))
	ereturn scalar me_y1_7  = -(`beta_min')*(normalden(`c2'-(`beta_min'*(0))-`xg_morex') - normalden(`c1'-(`beta_min'*(0))-`xg_morex'))
	ereturn scalar me_y1_8  = -(`beta_min')*(normalden(`c2'-(`beta_min'*(1))-`xg_morex') - normalden(`c1'-(`beta_min'*(1))-`xg_morex'))
	ereturn scalar me_y1_9  = -(`beta_min')*(normalden(`c2'-(`beta_min'*(2))-`xg_morex') - normalden(`c1'-(`beta_min'*(2))-`xg_morex'))
	ereturn scalar me_y1_10 = -(`beta_min')*(normalden(`c2'-(`beta_min'*(3))-`xg_morex') - normalden(`c1'-(`beta_min'*(3))-`xg_morex'))
	ereturn scalar me_y1_11 = -(`beta_min')*(normalden(`c2'-(`beta_min'*(4))-`xg_morex') - normalden(`c1'-(`beta_min'*(4))-`xg_morex'))
	
	ereturn scalar me_y2_1  = -(`beta_min')*(-normalden(`c2'-(`beta_min'*(-6))-`xg_morex'))
	ereturn scalar me_y2_2  = -(`beta_min')*(-normalden(`c2'-(`beta_min'*(-5))-`xg_morex'))
	ereturn scalar me_y2_3  = -(`beta_min')*(-normalden(`c2'-(`beta_min'*(-4))-`xg_morex'))
	ereturn scalar me_y2_4  = -(`beta_min')*(-normalden(`c2'-(`beta_min'*(-3))-`xg_morex'))
	ereturn scalar me_y2_5  = -(`beta_min')*(-normalden(`c2'-(`beta_min'*(-2))-`xg_morex'))
	ereturn scalar me_y2_6  = -(`beta_min')*(-normalden(`c2'-(`beta_min'*(-1))-`xg_morex'))
	ereturn scalar me_y2_7  = -(`beta_min')*(-normalden(`c2'-(`beta_min'*(0))-`xg_morex'))
	ereturn scalar me_y2_8  = -(`beta_min')*(-normalden(`c2'-(`beta_min'*(1))-`xg_morex'))
	ereturn scalar me_y2_9  = -(`beta_min')*(-normalden(`c2'-(`beta_min'*(2))-`xg_morex'))
	ereturn scalar me_y2_10 = -(`beta_min')*(-normalden(`c2'-(`beta_min'*(3))-`xg_morex'))
	ereturn scalar me_y2_11 = -(`beta_min')*(-normalden(`c2'-(`beta_min'*(4))-`xg_morex'))
	
end

bootstrap e(coef_min_cf), reps(100) nodrop: my_op_cf_min_x 
bootstrap e(coef_bal_cf), reps(100) nodrop: my_op_cf_min_x 
bootstrap e(coef_bos_cf), reps(100) nodrop: my_op_cf_min_x 
bootstrap e(coef_chi_cf), reps(100) nodrop: my_op_cf_min_x 
bootstrap e(coef_la_cf), reps(100) nodrop: my_op_cf_min_x 
bootstrap e(cut1), reps(100) nodrop: my_op_cf_min_x 
bootstrap e(cut2), reps(100) nodrop: my_op_cf_min_x 

*** conditional probabilities that Y = 0,1,2

forvalues y = 0(1)2{ 
	forvalues i = 1(1)11{
    
		scalar min_pr_y`y'_`i' = `i'-7

		bootstrap e(pr_y`y'_`i'), reps(100) nodrop: my_op_cf_min_x  
		scalar pr_y`y'_`i' 		= e(b)[1,1]
		scalar pr_y`y'_`i'_v	= e(V)[1,1]
	
		scalar pr_y`y'_`i'_lci = pr_y`y'_`i'-invnormal(0.975)*(sqrt(pr_y`y'_`i'_v))
		scalar pr_y`y'_`i'_uci = pr_y`y'_`i'+invnormal(0.975)*(sqrt(pr_y`y'_`i'_v))

	}

	matrix min_pr_y`y'_x = (scalar(min_pr_y`y'_1),scalar(pr_y`y'_1),scalar(pr_y`y'_1_lci),scalar(pr_y`y'_1_uci)\/*
				*/scalar(min_pr_y`y'_2),scalar(pr_y`y'_2),scalar(pr_y`y'_2_lci),scalar(pr_y`y'_2_uci)\/*
				*/scalar(min_pr_y`y'_3),scalar(pr_y`y'_3),scalar(pr_y`y'_3_lci),scalar(pr_y`y'_3_uci)\/*
				*/scalar(min_pr_y`y'_4),scalar(pr_y`y'_4),scalar(pr_y`y'_4_lci),scalar(pr_y`y'_4_uci)\/*
				*/scalar(min_pr_y`y'_5),scalar(pr_y`y'_5),scalar(pr_y`y'_5_lci),scalar(pr_y`y'_5_uci)\/*
				*/scalar(min_pr_y`y'_6),scalar(pr_y`y'_6),scalar(pr_y`y'_6_lci),scalar(pr_y`y'_6_uci)\/*
				*/scalar(min_pr_y`y'_7),scalar(pr_y`y'_7),scalar(pr_y`y'_7_lci),scalar(pr_y`y'_7_uci)\/*
				*/scalar(min_pr_y`y'_8),scalar(pr_y`y'_8),scalar(pr_y`y'_8_lci),scalar(pr_y`y'_8_uci)\/*
				*/scalar(min_pr_y`y'_9),scalar(pr_y`y'_9),scalar(pr_y`y'_9_lci),scalar(pr_y`y'_9_uci)\/*
				*/scalar(min_pr_y`y'_10),scalar(pr_y`y'_10),scalar(pr_y`y'_10_lci),scalar(pr_y`y'_10_uci)\/*
				*/scalar(min_pr_y`y'_11),scalar(pr_y`y'_11),scalar(pr_y`y'_11_lci),scalar(pr_y`y'_11_uci))

	matrix colnames min_pr_y`y'_x = minority_pr_y`y'_x min_probs_y`y'_x min_pr_lc_y`y'_x min_pr_uc_y`y'_x
	svmat min_pr_y`y'_x, names(col)
	twoway rcap min_pr_lc_y`y'_x min_pr_uc_y`y'_x minority_pr_y`y'_x, lcolor(black)||scatter min_probs_y`y'_x minority_pr_y`y'_x, mcolor(black)||line min_probs_y`y'_x minority_pr_y`y'_x, /*
	*/lcolor(black) title(" ") xtitle(" ") ytitle(" ") ylabel(-0.4(0.4)1.2) legend(off)
	
	graph export ${pathres}prob_CF_op_y`y'_x_min_med_ny_${date}.eps, replace
	
	scalar drop min_pr_y`y'_1 min_pr_y`y'_2 min_pr_y`y'_3 min_pr_y`y'_4 min_pr_y`y'_5 min_pr_y`y'_6 min_pr_y`y'_7 min_pr_y`y'_8 min_pr_y`y'_9 min_pr_y`y'_10 min_pr_y`y'_11 /*
			*/ pr_y`y'_1 pr_y`y'_2 pr_y`y'_3 pr_y`y'_4 pr_y`y'_5 pr_y`y'_6 pr_y`y'_7 pr_y`y'_8 pr_y`y'_9 pr_y`y'_10 pr_y`y'_11 /*
			*/ pr_y`y'_1_lci pr_y`y'_2_lci pr_y`y'_3_lci pr_y`y'_4_lci pr_y`y'_5_lci pr_y`y'_6_lci pr_y`y'_7_lci pr_y`y'_8_lci pr_y`y'_9_lci pr_y`y'_10_lci pr_y`y'_11_lci /*
			*/ pr_y`y'_1_uci pr_y`y'_2_uci pr_y`y'_3_uci pr_y`y'_4_uci pr_y`y'_5_uci pr_y`y'_6_uci pr_y`y'_7_uci pr_y`y'_8_uci pr_y`y'_9_uci pr_y`y'_10_uci pr_y`y'_11_uci  
}

*** conditional marginal effects for Y = 0,1,2

forvalues y = 0(1)2{ 
	forvalues i = 1(1)11{
    
		scalar min_me_y`y'_`i' = `i'-7

		bootstrap e(me_y`y'_`i'), reps(100) nodrop: my_op_cf_min_x 
		scalar me_y`y'_`i' 		= e(b)[1,1]
		scalar me_y`y'_`i'_v	= e(V)[1,1]
	
		scalar me_y`y'_`i'_lci = me_y`y'_`i'-invnormal(0.975)*(sqrt(me_y`y'_`i'_v))
		scalar me_y`y'_`i'_uci = me_y`y'_`i'+invnormal(0.975)*(sqrt(me_y`y'_`i'_v))

	}

	matrix min_me_y`y'_x = (scalar(min_me_y`y'_1),scalar(me_y`y'_1),scalar(me_y`y'_1_lci),scalar(me_y`y'_1_uci)\/*
				*/scalar(min_me_y`y'_2),scalar(me_y`y'_2),scalar(me_y`y'_2_lci),scalar(me_y`y'_2_uci)\/*
				*/scalar(min_me_y`y'_3),scalar(me_y`y'_3),scalar(me_y`y'_3_lci),scalar(me_y`y'_3_uci)\/*
				*/scalar(min_me_y`y'_4),scalar(me_y`y'_4),scalar(me_y`y'_4_lci),scalar(me_y`y'_4_uci)\/*
				*/scalar(min_me_y`y'_5),scalar(me_y`y'_5),scalar(me_y`y'_5_lci),scalar(me_y`y'_5_uci)\/*
				*/scalar(min_me_y`y'_6),scalar(me_y`y'_6),scalar(me_y`y'_6_lci),scalar(me_y`y'_6_uci)\/*
				*/scalar(min_me_y`y'_7),scalar(me_y`y'_7),scalar(me_y`y'_7_lci),scalar(me_y`y'_7_uci)\/*
				*/scalar(min_me_y`y'_8),scalar(me_y`y'_8),scalar(me_y`y'_8_lci),scalar(me_y`y'_8_uci)\/*
				*/scalar(min_me_y`y'_9),scalar(me_y`y'_9),scalar(me_y`y'_9_lci),scalar(me_y`y'_9_uci)\/*
				*/scalar(min_me_y`y'_10),scalar(me_y`y'_10),scalar(me_y`y'_10_lci),scalar(me_y`y'_10_uci)\/*
				*/scalar(min_me_y`y'_11),scalar(me_y`y'_11),scalar(me_y`y'_11_lci),scalar(me_y`y'_11_uci))

	matrix colnames min_me_y`y'_x = minority_me_y`y'_x min_meffects_y`y'_x min_me_lc_y`y'_x min_me_uc_y`y'_x
	svmat min_me_y`y'_x, names(col)
	twoway rcap min_me_lc_y`y'_x min_me_uc_y`y'_x minority_me_y`y'_x, lcolor(black)||scatter min_meffects_y`y'_x minority_me_y`y'_x, mcolor(black)||line min_meffects_y`y'_x minority_me_y`y'_x, lcolor(black) title(" ") xtitle(" ") ytitle(" ") ylabel(-0.2(0.1)0.2) legend(off)
	
	graph export ${pathres}me_CF_op_y`y'_x_min_med_ny_${date}.eps, replace
	
	scalar drop min_me_y`y'_1 min_me_y`y'_2 min_me_y`y'_3 min_me_y`y'_4 min_me_y`y'_5 min_me_y`y'_6 min_me_y`y'_7 min_me_y`y'_8 min_me_y`y'_9 min_me_y`y'_10 min_me_y`y'_11 /*
			*/ me_y`y'_1 me_y`y'_2 me_y`y'_3 me_y`y'_4 me_y`y'_5 me_y`y'_6 me_y`y'_7 me_y`y'_8 me_y`y'_9 me_y`y'_10 me_y`y'_11 /*
			*/ me_y`y'_1_lci me_y`y'_2_lci me_y`y'_3_lci me_y`y'_4_lci me_y`y'_5_lci me_y`y'_6_lci me_y`y'_7_lci me_y`y'_8_lci me_y`y'_9_lci me_y`y'_10_lci me_y`y'_11_lci /*
			*/ me_y`y'_1_uci me_y`y'_2_uci me_y`y'_3_uci me_y`y'_4_uci me_y`y'_5_uci me_y`y'_6_uci me_y`y'_7_uci me_y`y'_8_uci me_y`y'_9_uci me_y`y'_10_uci me_y`y'_11_uci 
}

** W = (hood poverty, hood minority), no covariates

capture noisily program drop my_op_cf_povmin_nox 

program my_op_cf_povmin_nox, eclass
	tempvar povmin_eq2 
	tempvar resid_pov_eq1 resid_min_eq1 resid_povmin_eq2 

	regress f_c9010t_perpov_dw_z i.ra_group#i.ra_site x_f_site_bal x_f_site_bos x_f_site_chi x_f_site_la [pw=$wt]
	predict `resid_pov_eq1', residuals 
	
	regress f_c9010t_pminorty_dw_z i.ra_group#i.ra_site x_f_site_bal x_f_site_bos x_f_site_chi x_f_site_la [pw=$wt]
	predict `resid_min_eq1', residuals 
	
	oprobit happy_scale012 `resid_pov_eq1' `resid_min_eq1' f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z x_f_site_bal x_f_site_bos x_f_site_chi x_f_site_la [pw=$wt]
	predict `povmin_eq2'
	gen `resid_povmin_eq2' = happy_scale012 - `povmin_eq2'
	
	correlate `resid_pov_eq1' `resid_povmin_eq2', covariance 
	local r_1 = r(cov_12)
	correlate `resid_min_eq1' `resid_povmin_eq2', covariance 
	local r_2 = r(cov_12)
	matrix r = (`r_1', `r_2')
		
	correlate `resid_pov_eq1' `resid_min_eq1', covariance
	matrix sigma = r(C)
	matrix mat_delta = r*inv(sigma)*r'
	
	local delta = (1 - el(mat_delta,1,1))^0.5
	
	ereturn scalar coef_povmin1_cf = _b[f_c9010t_perpov_dw_z]*`delta'
	ereturn scalar coef_povmin2_cf = _b[f_c9010t_pminorty_dw_z]*`delta'
	ereturn scalar coef_bal_cf = _b[x_f_site_bal]*`delta'
	ereturn scalar coef_bos_cf = _b[x_f_site_bos]*`delta'
	ereturn scalar coef_chi_cf = _b[x_f_site_chi]*`delta'
	ereturn scalar coef_la_cf  = _b[x_f_site_la]*`delta'
	ereturn scalar cut1		= _b[/:cut1]*`delta'
	ereturn scalar cut2		= _b[/:cut2]*`delta'

	local beta_povmin1  = _b[f_c9010t_perpov_dw_z]*`delta'
	local beta_povmin2  = _b[f_c9010t_pminorty_dw_z]*`delta'
	local gamma_bal = _b[x_f_site_bal]*`delta'
	local gamma_bos = _b[x_f_site_bos]*`delta'
	local gamma_chi = _b[x_f_site_chi]*`delta'
	local gamma_la  = _b[x_f_site_la]*`delta'
	
	local c1		= _b[/:cut1]*`delta'
	local c2		= _b[/:cut2]*`delta'
	
end

bootstrap e(coef_povmin1_cf), reps(100) nodrop: my_op_cf_povmin_nox 
bootstrap e(coef_povmin2_cf), reps(100) nodrop: my_op_cf_povmin_nox 
bootstrap e(coef_bal_cf), reps(100) nodrop: my_op_cf_povmin_nox 
bootstrap e(coef_bos_cf), reps(100) nodrop: my_op_cf_povmin_nox 
bootstrap e(coef_chi_cf), reps(100) nodrop: my_op_cf_povmin_nox 
bootstrap e(coef_la_cf), reps(100) nodrop: my_op_cf_povmin_nox 
bootstrap e(cut1), reps(100) nodrop: my_op_cf_povmin_nox 
bootstrap e(cut2), reps(100) nodrop: my_op_cf_povmin_nox 
	
** W = (hood poverty, hood minority), with covariates

capture noisily program drop my_op_cf_povmin_x 

program my_op_cf_povmin_x, eclass
	tempvar povmin_eq2 
	tempvar resid_pov_eq1 resid_min_eq1 resid_povmin_eq2 

	regress f_c9010t_perpov_dw_z i.ra_group#i.ra_site $site_covs x_f_release1 $x_covariates [pw=$wt] 
	predict `resid_pov_eq1', residuals 
	
	regress f_c9010t_pminorty_dw_z i.ra_group#i.ra_site $site_covs x_f_release1 $x_covariates [pw=$wt]
	predict `resid_min_eq1', residuals 
	
	oprobit happy_scale012 `resid_pov_eq1' `resid_min_eq1' f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z $site_covs x_f_release1 $x_covariates [pw=$wt] 
	predict `povmin_eq2'
	gen `resid_povmin_eq2' = happy_scale012 - `povmin_eq2'
	
	correlate `resid_pov_eq1' `resid_povmin_eq2', covariance 
	local r_1 = r(cov_12)
	correlate `resid_min_eq1' `resid_povmin_eq2', covariance 
	local r_2 = r(cov_12)
	matrix r = (`r_1', `r_2')
		
	correlate `resid_pov_eq1' `resid_min_eq1', covariance
	matrix sigma = r(C)
	matrix mat_delta = r*inv(sigma)*r'
	
	local delta = (1 - el(mat_delta,1,1))^0.5
	
	ereturn scalar coef_povmin1_cf = _b[f_c9010t_perpov_dw_z]*`delta'
	ereturn scalar coef_povmin2_cf = _b[f_c9010t_pminorty_dw_z]*`delta'
	ereturn scalar coef_bal_cf = _b[x_f_site_bal]*`delta'
	ereturn scalar coef_bos_cf = _b[x_f_site_bos]*`delta'
	ereturn scalar coef_chi_cf = _b[x_f_site_chi]*`delta'
	ereturn scalar coef_la_cf  = _b[x_f_site_la]*`delta'
	ereturn scalar cut1		= _b[/:cut1]*`delta'
	ereturn scalar cut2		= _b[/:cut2]*`delta'

	local beta_povmin1  = _b[f_c9010t_perpov_dw_z]*`delta'
	local beta_povmin2  = _b[f_c9010t_pminorty_dw_z]*`delta'
	
	local gamma_bal = _b[x_f_site_bal]*`delta'
	local gamma_bos = _b[x_f_site_bos]*`delta'
	local gamma_chi = _b[x_f_site_chi]*`delta'
	local gamma_la  = _b[x_f_site_la]*`delta'
	
	local c1		= _b[/:cut1]*`delta'
	local c2		= _b[/:cut2]*`delta'
	
	local g_x_f_release1 				= _b[x_f_release1]				*`delta'
    local g_tmp_female 					= _b[tmp_female]				*`delta'
    local g_x_rad_ad_le_35 				= _b[x_rad_ad_le_35]			*`delta'
    local g_x_rad_ad_36_40 				= _b[x_rad_ad_36_40]			*`delta'
    local g_x_rad_ad_41_45 				= _b[x_rad_ad_41_45]			*`delta'
    local g_x_rad_ad_46_50 				= _b[x_rad_ad_46_50]			*`delta'
	local g_x_rad_ad_ethrace_black_nh 	= _b[x_rad_ad_ethrace_black_nh]	*`delta'
    local g_x_rad_ad_ethrace_hisp 		= _b[x_rad_ad_ethrace_hisp]		*`delta'
    local g_x_f_ad_nevmarr 				= _b[x_f_ad_nevmarr]			*`delta'
    local g_x_f_ad_parentu18 			= _b[x_f_ad_parentu18]			*`delta'
    local g_x_f_ad_working 				= _b[x_f_ad_working]			*`delta'
    local g_x_f_ad_edinsch 				= _b[x_f_ad_edinsch]			*`delta'
    local g_x_f_ad_edgradhs 			= _b[x_f_ad_edgradhs]			*`delta'
    local g_x_f_ad_edged 				= _b[x_f_ad_edged]				*`delta'
    local g_x_f_hh_afdc 				= _b[x_f_hh_afdc]				*`delta'
	local g_rad_svy_bl_totincm_2009d 	= _b[rad_svy_bl_totincm_2009d]	*`delta'
    local g_x_f_hh_car 					= _b[x_f_hh_car]				*`delta'
    local g_x_f_hh_disabl 				= _b[x_f_hh_disabl]				*`delta'
    local g_x_f_hh_noteens 				= _b[x_f_hh_noteens]			*`delta'
    local g_x_f_hh_size2 				= _b[x_f_hh_size2]				*`delta'
    local g_x_f_hh_size3 				= _b[x_f_hh_size3]				*`delta'
    local g_x_f_hh_size4 				= _b[x_f_hh_size4]				*`delta'
    local g_x_f_hh_victim 				= _b[x_f_hh_victim]				*`delta'
    local g_x_f_hood_unsafenit 			= _b[x_f_hood_unsafenit]		*`delta'
    local g_x_f_hood_verydissat 		= _b[x_f_hood_verydissat]		*`delta'
    local g_x_f_hood_5y 				= _b[x_f_hood_5y]				*`delta'
    local g_x_f_hous_mov3tm 			= _b[x_f_hous_mov3tm]			*`delta'
    local g_x_f_hood_nofamily 			= _b[x_f_hood_nofamily]			*`delta'
    local g_x_f_hood_nofriend 			= _b[x_f_hood_nofriend]			*`delta'
    local g_x_f_hood_chat 				= _b[x_f_hood_chat]				*`delta'
    local g_x_f_hood_nbrkid 			= _b[x_f_hood_nbrkid]			*`delta'
    local g_x_f_hous_fndapt 			= _b[x_f_hous_fndapt]			*`delta'
    local g_x_f_hous_sec8bef 			= _b[x_f_hous_sec8bef]			*`delta'
    local g_x_f_hous_movdrgs 			= _b[x_f_hous_movdrgs]			*`delta'
    local g_x_f_hous_movschl 			= _b[x_f_hous_movschl]			*`delta'
    local g_cov_hous_movapt 			= _b[cov_hous_movapt]			*`delta'
    local g_cov_hous_movjob 			= _b[cov_hous_movjob]			*`delta'
	
	local xg_morex 	= (x_f_release1_med_x*`g_x_f_release1') + (tmp_female_med_x*`g_tmp_female') + (x_rad_ad_le_35_med_x*`g_x_rad_ad_le_35') + (x_rad_ad_36_40_med_x*`g_x_rad_ad_36_40') +(x_rad_ad_41_45_med_x*`g_x_rad_ad_41_45') + (x_rad_ad_46_50_med_x*`g_x_rad_ad_46_50') + (x_rad_ad_ethrace_black_nh_med_x*`g_x_rad_ad_ethrace_black_nh') + (x_rad_ad_ethrace_hisp_med_x*`g_x_rad_ad_ethrace_hisp') + (x_f_ad_nevmarr_med_x*`g_x_f_ad_nevmarr') + (x_f_ad_parentu18_med_x*`g_x_f_ad_parentu18') + (x_f_ad_working_med_x*`g_x_f_ad_working') + (x_f_ad_edinsch_med_x*`g_x_f_ad_edinsch') + (x_f_ad_edgradhs_med_x*`g_x_f_ad_edgradhs') + (x_f_ad_edged_med_x*`g_x_f_ad_edged') + (x_f_hh_afdc_med_x*`g_x_f_hh_afdc') + (rad_svy_bl_totincm_2009d_med_x*`g_rad_svy_bl_totincm_2009d') + (x_f_hh_car_med_x*`g_x_f_hh_car') + (x_f_hh_disabl_med_x*`g_x_f_hh_disabl') + (x_f_hh_noteens_med_x*`g_x_f_hh_noteens') + (x_f_hh_size2_med_x*`g_x_f_hh_size2') +(x_f_hh_size3_med_x*`g_x_f_hh_size3') + (x_f_hh_size4_med_x*`g_x_f_hh_size4') + (x_f_hh_victim_med_x*`g_x_f_hh_victim') + (x_f_hood_unsafenit_med_x*`g_x_f_hood_unsafenit') + (x_f_hood_verydissat_med_x*`g_x_f_hood_verydissat') + (x_f_hood_5y_med_x*`g_x_f_hood_5y') + (x_f_hous_mov3tm_med_x*`g_x_f_hous_mov3tm') + (x_f_hood_nofamily_med_x*`g_x_f_hood_nofamily') + (x_f_hood_nofriend_med_x*`g_x_f_hood_nofriend') + (x_f_hood_chat_med_x*`g_x_f_hood_chat') + (x_f_hood_nbrkid_med_x*`g_x_f_hood_nbrkid') + (x_f_hous_fndapt_med_x*`g_x_f_hous_fndapt') + (x_f_hous_sec8bef_med_x*`g_x_f_hous_sec8bef') + (x_f_hous_movdrgs_med_x*`g_x_f_hous_movdrgs') + (x_f_hous_movschl_med_x*`g_x_f_hous_movschl') + (cov_hous_movapt_med_x*`g_cov_hous_movapt') + (cov_hous_movjob_med_x*`g_cov_hous_movjob')
	
	ereturn scalar pr1_y0_1  = normal(`c1'-(`beta_povmin1'*(-6))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y0_2  = normal(`c1'-(`beta_povmin1'*(-5))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y0_3  = normal(`c1'-(`beta_povmin1'*(-4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y0_4  = normal(`c1'-(`beta_povmin1'*(-3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y0_5  = normal(`c1'-(`beta_povmin1'*(-2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y0_6  = normal(`c1'-(`beta_povmin1'*(-1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y0_7  = normal(`c1'-(`beta_povmin1'*(0))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y0_8  = normal(`c1'-(`beta_povmin1'*(1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y0_9  = normal(`c1'-(`beta_povmin1'*(2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y0_10 = normal(`c1'-(`beta_povmin1'*(3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y0_11 = normal(`c1'-(`beta_povmin1'*(4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	
	ereturn scalar pr1_y1_1  = normal(`c2'-(`beta_povmin1'*(-6))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin1'*(-6))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y1_2  = normal(`c2'-(`beta_povmin1'*(-5))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin1'*(-5))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y1_3  = normal(`c2'-(`beta_povmin1'*(-4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin1'*(-4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y1_4  = normal(`c2'-(`beta_povmin1'*(-3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin1'*(-3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y1_5  = normal(`c2'-(`beta_povmin1'*(-2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin1'*(-2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y1_6  = normal(`c2'-(`beta_povmin1'*(-1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin1'*(-1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y1_7  = normal(`c2'-(`beta_povmin1'*(0))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin1'*(0))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y1_8  = normal(`c2'-(`beta_povmin1'*(1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin1'*(1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y1_9  = normal(`c2'-(`beta_povmin1'*(2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin1'*(2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y1_10 = normal(`c2'-(`beta_povmin1'*(3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin1'*(3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y1_11 = normal(`c2'-(`beta_povmin1'*(4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin1'*(4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	
	ereturn scalar pr1_y2_1  = 1 - normal(`c2'-(`beta_povmin1'*(-6))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y2_2  = 1 - normal(`c2'-(`beta_povmin1'*(-5))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y2_3  = 1 - normal(`c2'-(`beta_povmin1'*(-4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y2_4  = 1 - normal(`c2'-(`beta_povmin1'*(-3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y2_5  = 1 - normal(`c2'-(`beta_povmin1'*(-2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y2_6  = 1 - normal(`c2'-(`beta_povmin1'*(-1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y2_7  = 1 - normal(`c2'-(`beta_povmin1'*(0))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y2_8  = 1 - normal(`c2'-(`beta_povmin1'*(1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y2_9  = 1 - normal(`c2'-(`beta_povmin1'*(2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y2_10 = 1 - normal(`c2'-(`beta_povmin1'*(3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr1_y2_11 = 1 - normal(`c2'-(`beta_povmin1'*(4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex')
	
	ereturn scalar pr2_y0_1  = normal(`c1'-(`beta_povmin2'*(-6))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y0_2  = normal(`c1'-(`beta_povmin2'*(-5))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y0_3  = normal(`c1'-(`beta_povmin2'*(-4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y0_4  = normal(`c1'-(`beta_povmin2'*(-3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y0_5  = normal(`c1'-(`beta_povmin2'*(-2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y0_6  = normal(`c1'-(`beta_povmin2'*(-1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y0_7  = normal(`c1'-(`beta_povmin2'*(0))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y0_8  = normal(`c1'-(`beta_povmin2'*(1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y0_9  = normal(`c1'-(`beta_povmin2'*(2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y0_10 = normal(`c1'-(`beta_povmin2'*(3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y0_11 = normal(`c1'-(`beta_povmin2'*(4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	
	ereturn scalar pr2_y1_1  = normal(`c2'-(`beta_povmin2'*(-6))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin2'*(-6))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y1_2  = normal(`c2'-(`beta_povmin2'*(-5))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin2'*(-5))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y1_3  = normal(`c2'-(`beta_povmin2'*(-4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin2'*(-4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y1_4  = normal(`c2'-(`beta_povmin2'*(-3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin2'*(-3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y1_5  = normal(`c2'-(`beta_povmin2'*(-2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin2'*(-2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y1_6  = normal(`c2'-(`beta_povmin2'*(-1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin2'*(-1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y1_7  = normal(`c2'-(`beta_povmin2'*(0))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin2'*(0))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y1_8  = normal(`c2'-(`beta_povmin2'*(1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin2'*(1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y1_9  = normal(`c2'-(`beta_povmin2'*(2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin2'*(2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y1_10 = normal(`c2'-(`beta_povmin2'*(3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin2'*(3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y1_11 = normal(`c2'-(`beta_povmin2'*(4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normal(`c1'-(`beta_povmin2'*(4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	
	ereturn scalar pr2_y2_1  = 1 - normal(`c2'-(`beta_povmin2'*(-6))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y2_2  = 1 - normal(`c2'-(`beta_povmin2'*(-5))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y2_3  = 1 - normal(`c2'-(`beta_povmin2'*(-4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y2_4  = 1 - normal(`c2'-(`beta_povmin2'*(-3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y2_5  = 1 - normal(`c2'-(`beta_povmin2'*(-2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y2_6  = 1 - normal(`c2'-(`beta_povmin2'*(-1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y2_7  = 1 - normal(`c2'-(`beta_povmin2'*(0))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y2_8  = 1 - normal(`c2'-(`beta_povmin2'*(1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y2_9  = 1 - normal(`c2'-(`beta_povmin2'*(2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y2_10 = 1 - normal(`c2'-(`beta_povmin2'*(3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	ereturn scalar pr2_y2_11 = 1 - normal(`c2'-(`beta_povmin2'*(4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex')
	
	ereturn scalar me1_y0_1  = -(`beta_povmin1')*(normalden(`c1'-(`beta_povmin1'*(-6))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y0_2  = -(`beta_povmin1')*(normalden(`c1'-(`beta_povmin1'*(-5))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y0_3  = -(`beta_povmin1')*(normalden(`c1'-(`beta_povmin1'*(-4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y0_4  = -(`beta_povmin1')*(normalden(`c1'-(`beta_povmin1'*(-3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y0_5  = -(`beta_povmin1')*(normalden(`c1'-(`beta_povmin1'*(-2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y0_6  = -(`beta_povmin1')*(normalden(`c1'-(`beta_povmin1'*(-1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y0_7  = -(`beta_povmin1')*(normalden(`c1'-(`beta_povmin1'*(0))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y0_8  = -(`beta_povmin1')*(normalden(`c1'-(`beta_povmin1'*(1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y0_9  = -(`beta_povmin1')*(normalden(`c1'-(`beta_povmin1'*(2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y0_10 = -(`beta_povmin1')*(normalden(`c1'-(`beta_povmin1'*(3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y0_11 = -(`beta_povmin1')*(normalden(`c1'-(`beta_povmin1'*(4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	
	ereturn scalar me1_y1_1  = -(`beta_povmin1')*(normalden(`c2'-(`beta_povmin1'*(-6))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin1'*(-6))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y1_2  = -(`beta_povmin1')*(normalden(`c2'-(`beta_povmin1'*(-5))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin1'*(-5))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y1_3  = -(`beta_povmin1')*(normalden(`c2'-(`beta_povmin1'*(-4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin1'*(-4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y1_4  = -(`beta_povmin1')*(normalden(`c2'-(`beta_povmin1'*(-3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin1'*(-3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y1_5  = -(`beta_povmin1')*(normalden(`c2'-(`beta_povmin1'*(-2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin1'*(-2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y1_6  = -(`beta_povmin1')*(normalden(`c2'-(`beta_povmin1'*(-1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin1'*(-1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y1_7  = -(`beta_povmin1')*(normalden(`c2'-(`beta_povmin1'*(0))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin1'*(0))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y1_8  = -(`beta_povmin1')*(normalden(`c2'-(`beta_povmin1'*(1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin1'*(1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y1_9  = -(`beta_povmin1')*(normalden(`c2'-(`beta_povmin1'*(2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin1'*(2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y1_10 = -(`beta_povmin1')*(normalden(`c2'-(`beta_povmin1'*(3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin1'*(3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y1_11 = -(`beta_povmin1')*(normalden(`c2'-(`beta_povmin1'*(4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin1'*(4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	
	ereturn scalar me1_y2_1  = -(`beta_povmin1')*(-normalden(`c2'-(`beta_povmin1'*(-6))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y2_2  = -(`beta_povmin1')*(-normalden(`c2'-(`beta_povmin1'*(-5))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y2_3  = -(`beta_povmin1')*(-normalden(`c2'-(`beta_povmin1'*(-4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y2_4  = -(`beta_povmin1')*(-normalden(`c2'-(`beta_povmin1'*(-3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y2_5  = -(`beta_povmin1')*(-normalden(`c2'-(`beta_povmin1'*(-2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y2_6  = -(`beta_povmin1')*(-normalden(`c2'-(`beta_povmin1'*(-1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y2_7  = -(`beta_povmin1')*(-normalden(`c2'-(`beta_povmin1'*(0))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y2_8  = -(`beta_povmin1')*(-normalden(`c2'-(`beta_povmin1'*(1))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y2_9  = -(`beta_povmin1')*(-normalden(`c2'-(`beta_povmin1'*(2))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y2_10 = -(`beta_povmin1')*(-normalden(`c2'-(`beta_povmin1'*(3))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me1_y2_11 = -(`beta_povmin1')*(-normalden(`c2'-(`beta_povmin1'*(4))-(`beta_povmin2'*f_c9010t_pminorty_dw_z_med_nox)-`xg_morex'))
	
	ereturn scalar me2_y0_1  = -(`beta_povmin2')*(normalden(`c1'-(`beta_povmin2'*(-6))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y0_2  = -(`beta_povmin2')*(normalden(`c1'-(`beta_povmin2'*(-5))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y0_3  = -(`beta_povmin2')*(normalden(`c1'-(`beta_povmin2'*(-4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y0_4  = -(`beta_povmin2')*(normalden(`c1'-(`beta_povmin2'*(-3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y0_5  = -(`beta_povmin2')*(normalden(`c1'-(`beta_povmin2'*(-2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y0_6  = -(`beta_povmin2')*(normalden(`c1'-(`beta_povmin2'*(-1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y0_7  = -(`beta_povmin2')*(normalden(`c1'-(`beta_povmin2'*(0))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y0_8  = -(`beta_povmin2')*(normalden(`c1'-(`beta_povmin2'*(1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y0_9  = -(`beta_povmin2')*(normalden(`c1'-(`beta_povmin2'*(2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y0_10 = -(`beta_povmin2')*(normalden(`c1'-(`beta_povmin2'*(3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y0_11 = -(`beta_povmin2')*(normalden(`c1'-(`beta_povmin2'*(4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	
	ereturn scalar me2_y1_1  = -(`beta_povmin2')*(normalden(`c2'-(`beta_povmin2'*(-6))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin2'*(-6))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y1_2  = -(`beta_povmin2')*(normalden(`c2'-(`beta_povmin2'*(-5))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin2'*(-5))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y1_3  = -(`beta_povmin2')*(normalden(`c2'-(`beta_povmin2'*(-4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin2'*(-4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y1_4  = -(`beta_povmin2')*(normalden(`c2'-(`beta_povmin2'*(-3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin2'*(-3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y1_5  = -(`beta_povmin2')*(normalden(`c2'-(`beta_povmin2'*(-2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin2'*(-2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y1_6  = -(`beta_povmin2')*(normalden(`c2'-(`beta_povmin2'*(-1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin2'*(-1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y1_7  = -(`beta_povmin2')*(normalden(`c2'-(`beta_povmin2'*(0))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin2'*(0))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y1_8  = -(`beta_povmin2')*(normalden(`c2'-(`beta_povmin2'*(1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin2'*(1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y1_9  = -(`beta_povmin2')*(normalden(`c2'-(`beta_povmin2'*(2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin2'*(2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y1_10 = -(`beta_povmin2')*(normalden(`c2'-(`beta_povmin2'*(3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin2'*(3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y1_11 = -(`beta_povmin2')*(normalden(`c2'-(`beta_povmin2'*(4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex') - /*
	*/ normalden(`c1'-(`beta_povmin2'*(4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	
	ereturn scalar me2_y2_1  = -(`beta_povmin2')*(-normalden(`c2'-(`beta_povmin2'*(-6))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y2_2  = -(`beta_povmin2')*(-normalden(`c2'-(`beta_povmin2'*(-5))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y2_3  = -(`beta_povmin2')*(-normalden(`c2'-(`beta_povmin2'*(-4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y2_4  = -(`beta_povmin2')*(-normalden(`c2'-(`beta_povmin2'*(-3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y2_5  = -(`beta_povmin2')*(-normalden(`c2'-(`beta_povmin2'*(-2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y2_6  = -(`beta_povmin2')*(-normalden(`c2'-(`beta_povmin2'*(-1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y2_7  = -(`beta_povmin2')*(-normalden(`c2'-(`beta_povmin2'*(0))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y2_8  = -(`beta_povmin2')*(-normalden(`c2'-(`beta_povmin2'*(1))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y2_9  = -(`beta_povmin2')*(-normalden(`c2'-(`beta_povmin2'*(2))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y2_10 = -(`beta_povmin2')*(-normalden(`c2'-(`beta_povmin2'*(3))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
	ereturn scalar me2_y2_11 = -(`beta_povmin2')*(-normalden(`c2'-(`beta_povmin2'*(4))-(`beta_povmin1'*f_c9010t_perpov_dw_z_med_nox)-`xg_morex'))
end

bootstrap e(coef_povmin1_cf), reps(100) nodrop: my_op_cf_povmin_x 
bootstrap e(coef_povmin2_cf), reps(100) nodrop: my_op_cf_povmin_x 
bootstrap e(coef_bal_cf), reps(100) nodrop: my_op_cf_povmin_x 
bootstrap e(coef_bos_cf), reps(100) nodrop: my_op_cf_povmin_x 
bootstrap e(coef_chi_cf), reps(100) nodrop: my_op_cf_povmin_x 
bootstrap e(coef_la_cf), reps(100) nodrop: my_op_cf_povmin_x 
bootstrap e(cut1), reps(100) nodrop: my_op_cf_povmin_x 
bootstrap e(cut2), reps(100) nodrop: my_op_cf_povmin_x 
	
*** conditional probabilities that Y = 0,1,2

forvalues y = 0(1)2{ 
	forvalues i = 1(1)11{
    
		scalar pov_pr_y`y'_`i' = `i'-7

		bootstrap e(pr1_y`y'_`i'), reps(100) nodrop: my_op_cf_povmin_x 
		scalar pr1_y`y'_`i' 		= e(b)[1,1]
		scalar pr1_y`y'_`i'_v	= e(V)[1,1]
	
		scalar pr1_y`y'_`i'_lci = pr1_y`y'_`i'-invnormal(0.975)*(sqrt(pr1_y`y'_`i'_v))
		scalar pr1_y`y'_`i'_uci = pr1_y`y'_`i'+invnormal(0.975)*(sqrt(pr1_y`y'_`i'_v))

	}

	matrix povmin1_pr_y`y'_x = (scalar(pov_pr_y`y'_1),scalar(pr1_y`y'_1),scalar(pr1_y`y'_1_lci),scalar(pr1_y`y'_1_uci)\/*
				*/scalar(pov_pr_y`y'_2),scalar(pr1_y`y'_2),scalar(pr1_y`y'_2_lci),scalar(pr1_y`y'_2_uci)\/*
				*/scalar(pov_pr_y`y'_3),scalar(pr1_y`y'_3),scalar(pr1_y`y'_3_lci),scalar(pr1_y`y'_3_uci)\/*
				*/scalar(pov_pr_y`y'_4),scalar(pr1_y`y'_4),scalar(pr1_y`y'_4_lci),scalar(pr1_y`y'_4_uci)\/*
				*/scalar(pov_pr_y`y'_5),scalar(pr1_y`y'_5),scalar(pr1_y`y'_5_lci),scalar(pr1_y`y'_5_uci)\/*
				*/scalar(pov_pr_y`y'_6),scalar(pr1_y`y'_6),scalar(pr1_y`y'_6_lci),scalar(pr1_y`y'_6_uci)\/*
				*/scalar(pov_pr_y`y'_7),scalar(pr1_y`y'_7),scalar(pr1_y`y'_7_lci),scalar(pr1_y`y'_7_uci)\/*
				*/scalar(pov_pr_y`y'_8),scalar(pr1_y`y'_8),scalar(pr1_y`y'_8_lci),scalar(pr1_y`y'_8_uci)\/*
				*/scalar(pov_pr_y`y'_9),scalar(pr1_y`y'_9),scalar(pr1_y`y'_9_lci),scalar(pr1_y`y'_9_uci)\/*
				*/scalar(pov_pr_y`y'_10),scalar(pr1_y`y'_10),scalar(pr1_y`y'_10_lci),scalar(pr1_y`y'_10_uci)\/*
				*/scalar(pov_pr_y`y'_11),scalar(pr1_y`y'_11),scalar(pr1_y`y'_11_lci),scalar(pr1_y`y'_11_uci))

	matrix colnames povmin1_pr_y`y'_x = poverty_pr1_y`y'_x povmin1_probs_y`y'_x povmin1_pr_lc_y`y'_x povmin1_pr_uc_y`y'_x
	svmat povmin1_pr_y`y'_x, names(col)
	twoway rcap povmin1_pr_lc_y`y'_x povmin1_pr_uc_y`y'_x poverty_pr1_y`y'_x, lcolor(black)||scatter povmin1_probs_y`y'_x poverty_pr1_y`y'_x, mcolor(black)||line povmin1_probs_y`y'_x poverty_pr1_y`y'_x, lcolor(black) title(" ") xtitle(" ") ytitle(" ") ylabel(-0.6(0.3)1.5) legend(off)
	
	graph export ${pathres}prob_CF_op_y`y'_x_povmin1_med_ny_${date}.eps, replace
	
	scalar drop pov_pr_y`y'_1 pov_pr_y`y'_2 pov_pr_y`y'_3 pov_pr_y`y'_4 pov_pr_y`y'_5 pov_pr_y`y'_6 pov_pr_y`y'_7 pov_pr_y`y'_8 pov_pr_y`y'_9 pov_pr_y`y'_10 pov_pr_y`y'_11 /*
			*/ pr1_y`y'_1 pr1_y`y'_2 pr1_y`y'_3 pr1_y`y'_4 pr1_y`y'_5 pr1_y`y'_6 pr1_y`y'_7 pr1_y`y'_8 pr1_y`y'_9 pr1_y`y'_10 pr1_y`y'_11 /*
			*/ pr1_y`y'_1_lci pr1_y`y'_2_lci pr1_y`y'_3_lci pr1_y`y'_4_lci pr1_y`y'_5_lci pr1_y`y'_6_lci pr1_y`y'_7_lci pr1_y`y'_8_lci pr1_y`y'_9_lci pr1_y`y'_10_lci pr1_y`y'_11_lci /*
			*/ pr1_y`y'_1_uci pr1_y`y'_2_uci pr1_y`y'_3_uci pr1_y`y'_4_uci pr1_y`y'_5_uci pr1_y`y'_6_uci pr1_y`y'_7_uci pr1_y`y'_8_uci pr1_y`y'_9_uci pr1_y`y'_10_uci pr1_y`y'_11_uci 
}

forvalues y = 0(1)2{ 
	forvalues i = 1(1)11{
    
		scalar min_pr_y`y'_`i' = `i'-7

		bootstrap e(pr2_y`y'_`i'), reps(100) nodrop: my_op_cf_povmin_x 
		scalar pr2_y`y'_`i' 		= e(b)[1,1]
		scalar pr2_y`y'_`i'_v	= e(V)[1,1]
	
		scalar pr2_y`y'_`i'_lci = pr2_y`y'_`i'-invnormal(0.975)*(sqrt(pr2_y`y'_`i'_v))
		scalar pr2_y`y'_`i'_uci = pr2_y`y'_`i'+invnormal(0.975)*(sqrt(pr2_y`y'_`i'_v))

	}

	matrix povmin2_pr_y`y'_x = (scalar(min_pr_y`y'_1),scalar(pr2_y`y'_1),scalar(pr2_y`y'_1_lci),scalar(pr2_y`y'_1_uci)\/*
				*/scalar(min_pr_y`y'_2),scalar(pr2_y`y'_2),scalar(pr2_y`y'_2_lci),scalar(pr2_y`y'_2_uci)\/*
				*/scalar(min_pr_y`y'_3),scalar(pr2_y`y'_3),scalar(pr2_y`y'_3_lci),scalar(pr2_y`y'_3_uci)\/*
				*/scalar(min_pr_y`y'_4),scalar(pr2_y`y'_4),scalar(pr2_y`y'_4_lci),scalar(pr2_y`y'_4_uci)\/*
				*/scalar(min_pr_y`y'_5),scalar(pr2_y`y'_5),scalar(pr2_y`y'_5_lci),scalar(pr2_y`y'_5_uci)\/*
				*/scalar(min_pr_y`y'_6),scalar(pr2_y`y'_6),scalar(pr2_y`y'_6_lci),scalar(pr2_y`y'_6_uci)\/*
				*/scalar(min_pr_y`y'_7),scalar(pr2_y`y'_7),scalar(pr2_y`y'_7_lci),scalar(pr2_y`y'_7_uci)\/*
				*/scalar(min_pr_y`y'_8),scalar(pr2_y`y'_8),scalar(pr2_y`y'_8_lci),scalar(pr2_y`y'_8_uci)\/*
				*/scalar(min_pr_y`y'_9),scalar(pr2_y`y'_9),scalar(pr2_y`y'_9_lci),scalar(pr2_y`y'_9_uci)\/*
				*/scalar(min_pr_y`y'_10),scalar(pr2_y`y'_10),scalar(pr2_y`y'_10_lci),scalar(pr2_y`y'_10_uci)\/*
				*/scalar(min_pr_y`y'_11),scalar(pr2_y`y'_11),scalar(pr2_y`y'_11_lci),scalar(pr2_y`y'_11_uci))

	matrix colnames povmin2_pr_y`y'_x = minority_pr2_y`y'_x povmin2_probs_y`y'_x povmin2_pr_lc_y`y'_x povmin2_pr_uc_y`y'_x
	svmat povmin2_pr_y`y'_x, names(col)
	twoway rcap povmin2_pr_lc_y`y'_x povmin2_pr_uc_y`y'_x minority_pr2_y`y'_x, lcolor(black)||scatter povmin2_probs_y`y'_x minority_pr2_y`y'_x, mcolor(black)||line povmin2_probs_y`y'_x minority_pr2_y`y'_x, lcolor(black) title(" ") xtitle(" ") ytitle(" ") ylabel(-0.6(0.3)1.5) legend(off)
	
	graph export ${pathres}prob_CF_op_y`y'_x_povmin2_med_ny_${date}.eps, replace
	
	scalar drop min_pr_y`y'_1 min_pr_y`y'_2 min_pr_y`y'_3 min_pr_y`y'_4 min_pr_y`y'_5 min_pr_y`y'_6 min_pr_y`y'_7 min_pr_y`y'_8 min_pr_y`y'_9 min_pr_y`y'_10 min_pr_y`y'_11 /*
			*/ pr2_y`y'_1 pr2_y`y'_2 pr2_y`y'_3 pr2_y`y'_4 pr2_y`y'_5 pr2_y`y'_6 pr2_y`y'_7 pr2_y`y'_8 pr2_y`y'_9 pr2_y`y'_10 pr2_y`y'_11 /*
			*/ pr2_y`y'_1_lci pr2_y`y'_2_lci pr2_y`y'_3_lci pr2_y`y'_4_lci pr2_y`y'_5_lci pr2_y`y'_6_lci pr2_y`y'_7_lci pr2_y`y'_8_lci pr2_y`y'_9_lci pr2_y`y'_10_lci pr2_y`y'_11_lci /*
			*/ pr2_y`y'_1_uci pr2_y`y'_2_uci pr2_y`y'_3_uci pr2_y`y'_4_uci pr2_y`y'_5_uci pr2_y`y'_6_uci pr2_y`y'_7_uci pr2_y`y'_8_uci pr2_y`y'_9_uci pr2_y`y'_10_uci pr2_y`y'_11_uci 
}

*** conditional marginal effects for Y = 0,1,2

forvalues y = 0(1)2{ 
	forvalues i = 1(1)11{
    
		scalar pov_me_y`y'_`i' = `i'-7

		bootstrap e(me1_y`y'_`i'), reps(100) nodrop: my_op_cf_povmin_x  
		scalar me1_y`y'_`i' 		= e(b)[1,1]
		scalar me1_y`y'_`i'_v	= e(V)[1,1]
	
		scalar me1_y`y'_`i'_lci = me1_y`y'_`i'-invnormal(0.975)*(sqrt(me1_y`y'_`i'_v))
		scalar me1_y`y'_`i'_uci = me1_y`y'_`i'+invnormal(0.975)*(sqrt(me1_y`y'_`i'_v))

	}

	matrix povmin1_me_y`y'_x = (scalar(pov_me_y`y'_1),scalar(me1_y`y'_1),scalar(me1_y`y'_1_lci),scalar(me1_y`y'_1_uci)\/*
				*/scalar(pov_me_y`y'_2),scalar(me1_y`y'_2),scalar(me1_y`y'_2_lci),scalar(me1_y`y'_2_uci)\/*
				*/scalar(pov_me_y`y'_3),scalar(me1_y`y'_3),scalar(me1_y`y'_3_lci),scalar(me1_y`y'_3_uci)\/*
				*/scalar(pov_me_y`y'_4),scalar(me1_y`y'_4),scalar(me1_y`y'_4_lci),scalar(me1_y`y'_4_uci)\/*
				*/scalar(pov_me_y`y'_5),scalar(me1_y`y'_5),scalar(me1_y`y'_5_lci),scalar(me1_y`y'_5_uci)\/*
				*/scalar(pov_me_y`y'_6),scalar(me1_y`y'_6),scalar(me1_y`y'_6_lci),scalar(me1_y`y'_6_uci)\/*
				*/scalar(pov_me_y`y'_7),scalar(me1_y`y'_7),scalar(me1_y`y'_7_lci),scalar(me1_y`y'_7_uci)\/*
				*/scalar(pov_me_y`y'_8),scalar(me1_y`y'_8),scalar(me1_y`y'_8_lci),scalar(me1_y`y'_8_uci)\/*
				*/scalar(pov_me_y`y'_9),scalar(me1_y`y'_9),scalar(me1_y`y'_9_lci),scalar(me1_y`y'_9_uci)\/*
				*/scalar(pov_me_y`y'_10),scalar(me1_y`y'_10),scalar(me1_y`y'_10_lci),scalar(me1_y`y'_10_uci)\/*
				*/scalar(pov_me_y`y'_11),scalar(me1_y`y'_11),scalar(me1_y`y'_11_lci),scalar(me1_y`y'_11_uci))

	matrix colnames povmin1_me_y`y'_x = poverty_me1_y`y'_x povmin1_meffects_y`y'_x povmin1_me_lc_y`y'_x povmin1_me_uc_y`y'_x
	svmat povmin1_me_y`y'_x, names(col)
	twoway rcap povmin1_me_lc_y`y'_x povmin1_me_uc_y`y'_x poverty_me1_y`y'_x, lcolor(black)||scatter povmin1_meffects_y`y'_x poverty_me1_y`y'_x, mcolor(black)||line povmin1_meffects_y`y'_x poverty_me1_y`y'_x, lcolor(black) title(" ") xtitle(" ") ytitle(" ") ylabel(-0.4(0.2)0.4) legend(off)
	
	graph export ${pathres}me_CF_op_y`y'_x_povmin1_med_ny_${date}.eps, replace
	
	scalar drop pov_me_y`y'_1 pov_me_y`y'_2 pov_me_y`y'_3 pov_me_y`y'_4 pov_me_y`y'_5 pov_me_y`y'_6 pov_me_y`y'_7 pov_me_y`y'_8 pov_me_y`y'_9 pov_me_y`y'_10 pov_me_y`y'_11 /*
			*/ me1_y`y'_1 me1_y`y'_2 me1_y`y'_3 me1_y`y'_4 me1_y`y'_5 me1_y`y'_6 me1_y`y'_7 me1_y`y'_8 me1_y`y'_9 me1_y`y'_10 me1_y`y'_11 /*
			*/ me1_y`y'_1_lci me1_y`y'_2_lci me1_y`y'_3_lci me1_y`y'_4_lci me1_y`y'_5_lci me1_y`y'_6_lci me1_y`y'_7_lci me1_y`y'_8_lci me1_y`y'_9_lci me1_y`y'_10_lci me1_y`y'_11_lci /*
			*/ me1_y`y'_1_uci me1_y`y'_2_uci me1_y`y'_3_uci me1_y`y'_4_uci me1_y`y'_5_uci me1_y`y'_6_uci me1_y`y'_7_uci me1_y`y'_8_uci me1_y`y'_9_uci me1_y`y'_10_uci me1_y`y'_11_uci 
}

forvalues y = 0(1)2{ 
	forvalues i = 1(1)11{
    
		scalar min_me_y`y'_`i' = `i'-7

		bootstrap e(me2_y`y'_`i'), reps(100) nodrop: my_op_cf_povmin_x  
		scalar me2_y`y'_`i' 		= e(b)[1,1]
		scalar me2_y`y'_`i'_v	= e(V)[1,1]
	
		scalar me2_y`y'_`i'_lci = me2_y`y'_`i'-invnormal(0.975)*(sqrt(me2_y`y'_`i'_v))
		scalar me2_y`y'_`i'_uci = me2_y`y'_`i'+invnormal(0.975)*(sqrt(me2_y`y'_`i'_v))

	}

	matrix povmin2_me_y`y'_x = (scalar(min_me_y`y'_1),scalar(me2_y`y'_1),scalar(me2_y`y'_1_lci),scalar(me2_y`y'_1_uci)\/*
				*/scalar(min_me_y`y'_2),scalar(me2_y`y'_2),scalar(me2_y`y'_2_lci),scalar(me2_y`y'_2_uci)\/*
				*/scalar(min_me_y`y'_3),scalar(me2_y`y'_3),scalar(me2_y`y'_3_lci),scalar(me2_y`y'_3_uci)\/*
				*/scalar(min_me_y`y'_4),scalar(me2_y`y'_4),scalar(me2_y`y'_4_lci),scalar(me2_y`y'_4_uci)\/*
				*/scalar(min_me_y`y'_5),scalar(me2_y`y'_5),scalar(me2_y`y'_5_lci),scalar(me2_y`y'_5_uci)\/*
				*/scalar(min_me_y`y'_6),scalar(me2_y`y'_6),scalar(me2_y`y'_6_lci),scalar(me2_y`y'_6_uci)\/*
				*/scalar(min_me_y`y'_7),scalar(me2_y`y'_7),scalar(me2_y`y'_7_lci),scalar(me2_y`y'_7_uci)\/*
				*/scalar(min_me_y`y'_8),scalar(me2_y`y'_8),scalar(me2_y`y'_8_lci),scalar(me2_y`y'_8_uci)\/*
				*/scalar(min_me_y`y'_9),scalar(me2_y`y'_9),scalar(me2_y`y'_9_lci),scalar(me2_y`y'_9_uci)\/*
				*/scalar(min_me_y`y'_10),scalar(me2_y`y'_10),scalar(me2_y`y'_10_lci),scalar(me2_y`y'_10_uci)\/*
				*/scalar(min_me_y`y'_11),scalar(me2_y`y'_11),scalar(me2_y`y'_11_lci),scalar(me2_y`y'_11_uci))

	matrix colnames povmin2_me_y`y'_x = minority_me2_y`y'_x povmin2_meffects_y`y'_x povmin2_me_lc_y`y'_x povmin2_me_uc_y`y'_x
	svmat povmin2_me_y`y'_x, names(col)
	twoway rcap povmin2_me_lc_y`y'_x povmin2_me_uc_y`y'_x minority_me2_y`y'_x, lcolor(black)||scatter povmin2_meffects_y`y'_x minority_me2_y`y'_x, mcolor(black)||line povmin2_meffects_y`y'_x minority_me2_y`y'_x, lcolor(black) title(" ") xtitle(" ") ytitle(" ") ylabel(-0.4(0.2)0.4) legend(off)
	
	graph export ${pathres}me_CF_op_y`y'_x_povmin2_med_ny_${date}.eps, replace
	
	scalar drop min_me_y`y'_1 min_me_y`y'_2 min_me_y`y'_3 min_me_y`y'_4 min_me_y`y'_5 min_me_y`y'_6 min_me_y`y'_7 min_me_y`y'_8 min_me_y`y'_9 min_me_y`y'_10 min_me_y`y'_11 /*
			*/ me2_y`y'_1 me2_y`y'_2 me2_y`y'_3 me2_y`y'_4 me2_y`y'_5 me2_y`y'_6 me2_y`y'_7 me2_y`y'_8 me2_y`y'_9 me2_y`y'_10 me2_y`y'_11 /*
			*/ me2_y`y'_1_lci me2_y`y'_2_lci me2_y`y'_3_lci me2_y`y'_4_lci me2_y`y'_5_lci me2_y`y'_6_lci me2_y`y'_7_lci me2_y`y'_8_lci me2_y`y'_9_lci me2_y`y'_10_lci me2_y`y'_11_lci /*
			*/ me2_y`y'_1_uci me2_y`y'_2_uci me2_y`y'_3_uci me2_y`y'_4_uci me2_y`y'_5_uci me2_y`y'_6_uci me2_y`y'_7_uci me2_y`y'_8_uci me2_y`y'_9_uci me2_y`y'_10_uci me2_y`y'_11_uci 
}