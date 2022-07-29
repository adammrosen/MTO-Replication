
* This file carries out IM tests after estimating MLE (CT model) 

* start with the case in which there is a single endogenous variable

* get estimates of a CT model in which the endogenous variable is neighborhood poverty

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs f_c9010t_perpov_dw_z, nocons) (cut1:) (lndiff:) (xb2:f_c9010t_perpov_dw_z = interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml maximize
est store pov_t_nox

* 1st and 2nd derivatives 

scalar infinity = 1000000000

gen V = f_c9010t_perpov_dw_z - ([xb2]_cons+/*
		*/(interaction_1_1*[xb2]interaction_1_1)+(interaction_1_2*[xb2]interaction_1_2)+(interaction_1_3*[xb2]interaction_1_3)+(interaction_1_4*[xb2]interaction_1_4)+/*
		*/(interaction_1_5*[xb2]interaction_1_5)+(interaction_2_1*[xb2]interaction_2_1)+(interaction_2_2*[xb2]interaction_2_2)+(interaction_2_3*[xb2]interaction_2_3)+/*
		*/(interaction_2_4*[xb2]interaction_2_4)+(interaction_2_5*[xb2]interaction_2_5)+(interaction_3_1*[xb2]interaction_3_1)+(interaction_3_2*[xb2]interaction_3_2)+/*
		*/(interaction_3_3*[xb2]interaction_3_3)+(interaction_3_4*[xb2]interaction_3_4))

gen op0 = -infinity
gen op1 = ([cut1]_cons-f_c9010t_perpov_dw_z*[xb1]f_c9010t_perpov_dw_z-x_f_site_balt*[xb1]x_f_site_balt-x_f_site_bos*[xb1]x_f_site_bos-x_f_site_chi*[xb1]x_f_site_chi-/*
																		*/x_f_site_la*[xb1]x_f_site_la-(([rho]_cons/exp([lnsigma]_cons))*V))
gen op2 = (exp([lndiff]_cons)+[cut1]_cons-f_c9010t_perpov_dw_z*[xb1]f_c9010t_perpov_dw_z-x_f_site_balt*[xb1]x_f_site_balt-x_f_site_bos*[xb1]x_f_site_bos-/*
																						*/x_f_site_chi*[xb1]x_f_site_chi-x_f_site_la*[xb1]x_f_site_la-(([rho]_cons/exp([lnsigma]_cons))*V))
gen op3 = infinity

scalar  d_rho = ((1-([rho]_cons^2))^0.5)

gen 	big_phi = normal(op1/d_rho)- normal(op0/d_rho) 	if happy_scale012 == 0									
replace big_phi = normal(op2/d_rho)- normal(op1/d_rho) 	if happy_scale012 == 1
replace big_phi = normal(op3/d_rho)- normal(op2/d_rho) 	if happy_scale012 == 2

gen 	small_phi = normalden(op1/d_rho)- normalden(op0/d_rho)	if happy_scale012 == 0									
replace small_phi = normalden(op2/d_rho)- normalden(op1/d_rho) 	if happy_scale012 == 1
replace small_phi = normalden(op3/d_rho)- normalden(op2/d_rho) 	if happy_scale012 == 2

gen f_beta 		= (-f_c9010t_perpov_dw_z/d_rho)*(small_phi/big_phi) 	

gen f_balt 		= (-x_f_site_balt/d_rho)*(small_phi/big_phi)
gen f_bos 		= (-x_f_site_bos/d_rho)*(small_phi/big_phi)
gen f_chi 		= (-x_f_site_chi/d_rho)*(small_phi/big_phi)
gen f_la 		= (-x_f_site_la/d_rho)*(small_phi/big_phi)

gen 	f_c2 	= 0 										if happy_scale012 == 0	
replace f_c2 	= (1/d_rho)*(normalden(op2/d_rho)/big_phi) 	if happy_scale012 == 1	
replace f_c2 	= (-1/d_rho)*(normalden(op2/d_rho)/big_phi) if happy_scale012 == 2	

gen 	f_c1 	= (1/d_rho)*(normalden(op1/d_rho)/big_phi) 	if happy_scale012 == 0	
replace f_c1 	= (-1/d_rho)*(normalden(op1/d_rho)/big_phi) if happy_scale012 == 1	
replace f_c1 	= 0 										if happy_scale012 == 2	

gen X_op0 = (-V/(exp([lnsigma]_cons)*d_rho))+(([rho]_cons/(d_rho^3))*op0)
gen X_op1 = (-V/(exp([lnsigma]_cons)*d_rho))+(([rho]_cons/(d_rho^3))*op1)
gen X_op2 = (-V/(exp([lnsigma]_cons)*d_rho))+(([rho]_cons/(d_rho^3))*op2)
gen X_op3 = (-V/(exp([lnsigma]_cons)*d_rho))+(([rho]_cons/(d_rho^3))*op3)

gen f_rho 		= ((normalden(op1/d_rho)*X_op1)- ((normalden(op0/d_rho)*X_op0)))/big_phi if happy_scale012 == 0	
replace f_rho 	= ((normalden(op2/d_rho)*X_op2)- ((normalden(op1/d_rho)*X_op1)))/big_phi if happy_scale012 == 1	
replace f_rho 	= ((normalden(op3/d_rho)*X_op3)- ((normalden(op2/d_rho)*X_op2)))/big_phi if happy_scale012 == 2	

gen f_sigma = ((([rho]_cons*V)/((exp([lnsigma]_cons)^2)*d_rho))*(small_phi/big_phi))+((1/exp([lnsigma]_cons))*(((V^2)/(exp([lnsigma]_cons)^2))-1))

gen x_f_site_ny = (ra_site==5)

gen f_deltaZ_11 = ((([rho]_cons*(x_f_site_balt*ra_grp_exp))/((exp([lnsigma]_cons)*d_rho)))*(small_phi/big_phi))+((x_f_site_balt*ra_grp_exp)*(V/(exp([lnsigma]_cons)^2)))
gen f_deltaZ_12 = ((([rho]_cons*(x_f_site_bos*ra_grp_exp))/((exp([lnsigma]_cons)*d_rho)))*(small_phi/big_phi))+((x_f_site_bos*ra_grp_exp)*(V/(exp([lnsigma]_cons)^2)))
gen f_deltaZ_13 = ((([rho]_cons*(x_f_site_chi*ra_grp_exp))/((exp([lnsigma]_cons)*d_rho)))*(small_phi/big_phi))+((x_f_site_chi*ra_grp_exp)*(V/(exp([lnsigma]_cons)^2)))
gen f_deltaZ_14 = ((([rho]_cons*(x_f_site_la*ra_grp_exp))/((exp([lnsigma]_cons)*d_rho)))*(small_phi/big_phi))+((x_f_site_la*ra_grp_exp)*(V/(exp([lnsigma]_cons)^2)))
gen f_deltaZ_15 = ((([rho]_cons*(x_f_site_ny*ra_grp_exp))/((exp([lnsigma]_cons)*d_rho)))*(small_phi/big_phi))+((x_f_site_ny*ra_grp_exp)*(V/(exp([lnsigma]_cons)^2)))

gen f_deltaZ_21 = ((([rho]_cons*(x_f_site_balt*ra_grp_s8))/((exp([lnsigma]_cons)*d_rho)))*(small_phi/big_phi))+((x_f_site_balt*ra_grp_s8)*(V/(exp([lnsigma]_cons)^2)))
gen f_deltaZ_22 = ((([rho]_cons*(x_f_site_bos*ra_grp_s8))/((exp([lnsigma]_cons)*d_rho)))*(small_phi/big_phi))+((x_f_site_bos*ra_grp_s8)*(V/(exp([lnsigma]_cons)^2)))
gen f_deltaZ_23 = ((([rho]_cons*(x_f_site_chi*ra_grp_s8))/((exp([lnsigma]_cons)*d_rho)))*(small_phi/big_phi))+((x_f_site_chi*ra_grp_s8)*(V/(exp([lnsigma]_cons)^2)))
gen f_deltaZ_24 = ((([rho]_cons*(x_f_site_la*ra_grp_s8))/((exp([lnsigma]_cons)*d_rho)))*(small_phi/big_phi))+((x_f_site_la*ra_grp_s8)*(V/(exp([lnsigma]_cons)^2)))
gen f_deltaZ_25 = ((([rho]_cons*(x_f_site_ny*ra_grp_s8))/((exp([lnsigma]_cons)*d_rho)))*(small_phi/big_phi))+((x_f_site_ny*ra_grp_s8)*(V/(exp([lnsigma]_cons)^2)))

gen f_deltaZ_31 = ((([rho]_cons*(x_f_site_balt*ra_grp_control))/((exp([lnsigma]_cons)*d_rho)))*(small_phi/big_phi))+((x_f_site_balt*ra_grp_control)*(V/(exp([lnsigma]_cons)^2)))
gen f_deltaZ_32 = ((([rho]_cons*(x_f_site_bos*ra_grp_control))/((exp([lnsigma]_cons)*d_rho)))*(small_phi/big_phi))+((x_f_site_bos*ra_grp_control)*(V/(exp([lnsigma]_cons)^2)))
gen f_deltaZ_33 = ((([rho]_cons*(x_f_site_chi*ra_grp_control))/((exp([lnsigma]_cons)*d_rho)))*(small_phi/big_phi))+((x_f_site_chi*ra_grp_control)*(V/(exp([lnsigma]_cons)^2)))
gen f_deltaZ_34 = ((([rho]_cons*(x_f_site_la*ra_grp_control))/((exp([lnsigma]_cons)*d_rho)))*(small_phi/big_phi))+((x_f_site_la*ra_grp_control)*(V/(exp([lnsigma]_cons)^2)))

gen f_alpha 	= (([rho]_cons)/(exp([lnsigma]_cons)*d_rho))*(small_phi/big_phi)+((V/(exp([lnsigma]_cons)^2)))

gen N = (-f_c9010t_perpov_dw_z/d_rho)*small_phi 

gen diff_N_inner 		= ((normalden(op1/d_rho))*(op1/d_rho)) - ((normalden(op0/d_rho))*(op0/d_rho)) if happy_scale012 == 0									
replace diff_N_inner 	= ((normalden(op2/d_rho))*(op2/d_rho)) - ((normalden(op1/d_rho))*(op1/d_rho)) if happy_scale012 == 1
replace diff_N_inner 	= ((normalden(op3/d_rho))*(op3/d_rho)) - ((normalden(op2/d_rho))*(op2/d_rho)) if happy_scale012 == 2

gen diff_N_beta = (-(f_c9010t_perpov_dw_z^2)/(d_rho^2))*(diff_N_inner) 

gen f_beta_beta = ((diff_N_beta*big_phi)-(N^2))/(big_phi^2)

gen A_balt = -(x_f_site_balt/(d_rho))*(small_phi)
gen diff_A_balt = (-(x_f_site_balt^2)/(d_rho^2))*(diff_N_inner) 
gen f_balt_balt = ((diff_A_balt*big_phi)-(A_balt^2))/(big_phi^2)

gen A_bos = -(x_f_site_bos/(d_rho))*(small_phi)
gen diff_A_bos = (-(x_f_site_bos^2)/(d_rho^2))*(diff_N_inner) 
gen f_bos_bos = ((diff_A_bos*big_phi)-(A_bos^2))/(big_phi^2)

gen A_chi = -(x_f_site_chi/(d_rho))*(small_phi)
gen diff_A_chi = (-(x_f_site_chi^2)/(d_rho^2))*(diff_N_inner) 
gen f_chi_chi = ((diff_A_chi*big_phi)-(A_chi^2))/(big_phi^2)

gen A_la = -(x_f_site_la/(d_rho))*(small_phi)
gen diff_A_la = (-(x_f_site_la^2)/(d_rho^2))*(diff_N_inner) 
gen f_la_la = ((diff_A_la*big_phi)-(A_la^2))/(big_phi^2)

gen diff_N_balt 	= (-(f_c9010t_perpov_dw_z*x_f_site_balt)/(d_rho^2))*(diff_N_inner) 
gen diff_N_bos 		= (-(f_c9010t_perpov_dw_z*x_f_site_bos)/(d_rho^2))*(diff_N_inner) 
gen diff_N_chi 		= (-(f_c9010t_perpov_dw_z*x_f_site_chi)/(d_rho^2))*(diff_N_inner) 
gen diff_N_la 		= (-(f_c9010t_perpov_dw_z*x_f_site_la)/(d_rho^2))*(diff_N_inner) 

gen diff_D_balt 	= -(x_f_site_balt/(d_rho))*(small_phi)
gen diff_D_bos 		= -(x_f_site_bos/(d_rho))*(small_phi)
gen diff_D_chi 		= -(x_f_site_chi/(d_rho))*(small_phi)
gen diff_D_la 		= -(x_f_site_la/(d_rho))*(small_phi)

gen f_balt_beta 	= ((diff_N_balt*big_phi)-(N*diff_D_balt))/(big_phi^2)
gen f_bos_beta 		= ((diff_N_bos*big_phi)-(N*diff_D_bos))/(big_phi^2)
gen f_chi_beta 		= ((diff_N_chi*big_phi)-(N*diff_D_chi))/(big_phi^2)
gen f_la_beta 		= ((diff_N_la*big_phi)-(N*diff_D_la))/(big_phi^2)

gen combine_beta_beta 	= (f_beta*f_beta) + f_beta_beta 
gen combine_balt_balt 	= (f_balt*f_balt) + f_balt_balt 
gen combine_bos_bos 	= (f_bos*f_bos) + f_bos_bos 
gen combine_chi_chi 	= (f_chi*f_chi) + f_chi_chi 
gen combine_la_la 		= (f_la*f_la) + f_la_la 

gen combine_beta_balt 	= (f_beta*f_balt) + f_balt_beta 
gen combine_beta_bos 	= (f_beta*f_bos) + f_bos_beta 
gen combine_beta_chi 	= (f_beta*f_chi) + f_chi_beta 
gen combine_beta_la 	= (f_beta*f_la) + f_la_beta 

* now run regressions of 1st and 2nd derivatives on column of ones to get the IM test statistic  

gen ones = 1

regress ones f_beta f_balt f_bos f_chi f_la f_c2 f_c1 f_alpha f_deltaZ_11 f_deltaZ_12 f_deltaZ_13 f_deltaZ_14 f_deltaZ_15 /*
			*/f_deltaZ_21 f_deltaZ_22 f_deltaZ_23 f_deltaZ_24 f_deltaZ_25 f_deltaZ_31 f_deltaZ_32 f_deltaZ_33 f_deltaZ_34 /*
			*/f_rho f_sigma, noconstant

* n=3263 and R-square=0.0023 in the regression below which gives an IM test statistic of 7.5049 or approximately 7.5. 
* The p-value reported in the paper comes from a chi-square distribution with 3 dof and equals 0.0576 

regress ones f_beta f_balt f_bos f_chi f_la f_c2 f_c1 f_alpha f_deltaZ_11 f_deltaZ_12 f_deltaZ_13 f_deltaZ_14 f_deltaZ_15 /*
			*/f_deltaZ_21 f_deltaZ_22 f_deltaZ_23 f_deltaZ_24 f_deltaZ_25 f_deltaZ_31 f_deltaZ_32 f_deltaZ_33 f_deltaZ_34 /*
			*/f_rho f_sigma combine_beta_beta combine_balt_balt combine_beta_balt, noconstant

* n=3263 and R-square=0.0051 in the regression below which gives an IM test statistic of 16.6413 or approximately 16.6. 
* The p-value reported in the paper comes from a chi-square distribution with 9 dof and equals 0.0554 
			
regress ones f_beta f_balt f_bos f_chi f_la f_c2 f_c1 f_alpha f_deltaZ_11 f_deltaZ_12 f_deltaZ_13 f_deltaZ_14 f_deltaZ_15 /*
			*/f_deltaZ_21 f_deltaZ_22 f_deltaZ_23 f_deltaZ_24 f_deltaZ_25 f_deltaZ_31 f_deltaZ_32 f_deltaZ_33 f_deltaZ_34 /*
			*/f_rho f_sigma combine_beta_beta combine_balt_balt combine_bos_bos combine_chi_chi combine_la_la combine_beta_balt /*
			*/combine_beta_bos combine_beta_chi combine_beta_la, noconstant

* second case in which there are multiple endogenous variables

* get estimates of a CT model in which the endogenous variables are neighborhood poverty & minority

ml model lf myoprobit3_lf (xb1:happy_scale012 = $site_covs f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z, nocons) (cut1:) (lndiff:) /*
						*/(xb2:f_c9010t_perpov_dw_z = interaction*) (lnvar_v1:) (rho_uv1:) /*
						*/(xb3:f_c9010t_pminorty_dw_z = interaction*) (lnvar_v2:) (rho_uv2:) (cov_v1v2:) [pw=$wt], /*
						*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
						*/diparm(lnvar_v1, exp label(var_v1)) /*
						*/diparm(lnvar_v2, exp label(var_v2)) 
ml init xb1:x_f_site_bal=0.3928 x_f_site_bos=0.4820 x_f_site_chi=0.2571 x_f_site_la=0.0768 f_c9010t_perpov_dw_z=-0.2857 f_c9010t_pminorty_dw_z=0.3206 cut1:_cons=-0.4024
ml maximize, difficult

* 1st and 2nd derivatives 

gen V1 	= f_c9010t_perpov_dw_z - ([xb2]_cons+(interaction_1_1*[xb2]interaction_1_1)+(interaction_1_2*[xb2]interaction_1_2)+(interaction_1_3*[xb2]interaction_1_3)+/*
											*/(interaction_1_4*[xb2]interaction_1_4)+(interaction_1_5*[xb2]interaction_1_5)+(interaction_2_1*[xb2]interaction_2_1)+/*
											*/(interaction_2_2*[xb2]interaction_2_2)+(interaction_2_3*[xb2]interaction_2_3)+(interaction_2_4*[xb2]interaction_2_4)+/*
											*/(interaction_2_5*[xb2]interaction_2_5)+(interaction_3_1*[xb2]interaction_3_1)+(interaction_3_2*[xb2]interaction_3_2)+/*
											*/(interaction_3_3*[xb2]interaction_3_3)+(interaction_3_4*[xb2]interaction_3_4))

gen V2 	= f_c9010t_pminorty_dw_z - ([xb3]_cons+(interaction_1_1*[xb3]interaction_1_1)+(interaction_1_2*[xb3]interaction_1_2)+(interaction_1_3*[xb3]interaction_1_3)+/*
											*/(interaction_1_4*[xb3]interaction_1_4)+(interaction_1_5*[xb3]interaction_1_5)+(interaction_2_1*[xb3]interaction_2_1)+/*
											*/(interaction_2_2*[xb3]interaction_2_2)+(interaction_2_3*[xb3]interaction_2_3)+(interaction_2_4*[xb3]interaction_2_4)+/*
											*/(interaction_2_5*[xb3]interaction_2_5)+(interaction_3_1*[xb3]interaction_3_1)+(interaction_3_2*[xb3]interaction_3_2)+/*
											*/(interaction_3_3*[xb3]interaction_3_3)+(interaction_3_4*[xb3]interaction_3_4))

scalar 	det = (exp([lnvar_v1]_cons)*exp([lnvar_v2]_cons))-(([cov_v1v2]_cons)^2) 

scalar 	sigma_v 	= ((det-(([rho_uv1]_cons^2)*exp([lnvar_v2]_cons))+(((2*[rho_uv1]_cons)*[rho_uv2]_cons)*[cov_v1v2]_cons)-(([rho_uv2]_cons^2)*exp([lnvar_v1]_cons)))/det)^0.5 
gen 	mu_v_num 	= (V1*(([rho_uv1]_cons*exp([lnvar_v2]_cons))-([rho_uv2]_cons*[cov_v1v2]_cons))+V2*(([rho_uv2]_cons*exp([lnvar_v1]_cons))-([rho_uv1]_cons*[cov_v1v2]_cons))) 
gen 	mu_v 		= mu_v_num/det 

replace op0 = -infinity
replace op1 = [cut1]_cons-(f_c9010t_perpov_dw_z*[xb1]f_c9010t_perpov_dw_z)-(f_c9010t_pminorty_dw_z*[xb1]f_c9010t_pminorty_dw_z)-/*
						*/(x_f_site_balt*[xb1]x_f_site_balt)-(x_f_site_bos*[xb1]x_f_site_bos)-(x_f_site_chi*[xb1]x_f_site_chi)-(x_f_site_la*[xb1]x_f_site_la)-mu_v
replace op2 = exp([lndiff]_cons)+[cut1]_cons-(f_c9010t_perpov_dw_z*[xb1]f_c9010t_perpov_dw_z)-(f_c9010t_pminorty_dw_z*[xb1]f_c9010t_pminorty_dw_z)-/*
						*/(x_f_site_balt*[xb1]x_f_site_balt)-(x_f_site_bos*[xb1]x_f_site_bos)-(x_f_site_chi*[xb1]x_f_site_chi)-(x_f_site_la*[xb1]x_f_site_la)-mu_v
replace op3 = infinity

replace big_phi 	= normal(op1/sigma_v)   - normal(op0/sigma_v) 	if happy_scale012 == 0									
replace big_phi 	= normal(op2/sigma_v)   - normal(op1/sigma_v) 	if happy_scale012 == 1
replace big_phi 	= normal(op3/sigma_v)   - normal(op2/sigma_v) 	if happy_scale012 == 2

replace small_phi 	= normalden(op1/sigma_v)- normalden(op0/sigma_v)	if happy_scale012 == 0									
replace small_phi 	= normalden(op2/sigma_v)- normalden(op1/sigma_v) 	if happy_scale012 == 1
replace small_phi 	= normalden(op3/sigma_v)- normalden(op2/sigma_v) 	if happy_scale012 == 2

gen fm_beta1 = (-f_c9010t_perpov_dw_z/sigma_v)*(small_phi/big_phi) 	
gen fm_beta2 = (-f_c9010t_pminorty_dw_z/sigma_v)*(small_phi/big_phi) 	

gen fm_balt  = (-x_f_site_balt/sigma_v)*(small_phi/big_phi)
gen fm_bos 	 = (-x_f_site_bos/sigma_v)*(small_phi/big_phi)
gen fm_chi 	 = (-x_f_site_chi/sigma_v)*(small_phi/big_phi)
gen fm_la 	 = (-x_f_site_la/sigma_v)*(small_phi/big_phi)

gen 	fm_c2 = 0 													if happy_scale012 == 0	
replace fm_c2 = (1/sigma_v)*(normalden(op2/sigma_v)/big_phi) 		if happy_scale012 == 1	
replace fm_c2 = (-1/sigma_v)*(normalden(op2/sigma_v)/big_phi) 		if happy_scale012 == 2	

gen 	fm_c1 = (1/sigma_v)*(normalden(op1/sigma_v)/big_phi) 		if happy_scale012 == 0	
replace fm_c1 = (-1/sigma_v)*(normalden(op1/sigma_v)/big_phi) 		if happy_scale012 == 1	
replace fm_c1 = 0 													if happy_scale012 == 2	

scalar 	fm_V1_num			= ((([rho_uv1]_cons*exp([lnvar_v2]_cons))-([rho_uv2]_cons*[cov_v1v2]_cons))/det)
gen 	fm_W_condXZ_diffV1	= ((V1*exp([lnvar_v2]_cons))-(V2*[cov_v1v2]_cons))/det

gen 	fm_alpha1 	  = (((1/sigma_v)*(fm_V1_num))*(small_phi/big_phi))+fm_W_condXZ_diffV1

gen 	fm_deltaZ1_11 = ((((x_f_site_balt*ra_grp_exp)/sigma_v)*fm_V1_num)*(small_phi/big_phi))+((x_f_site_balt*ra_grp_exp)*fm_W_condXZ_diffV1)
gen 	fm_deltaZ1_12 = ((((x_f_site_bos*ra_grp_exp)/sigma_v)*fm_V1_num)*(small_phi/big_phi)) +((x_f_site_bos*ra_grp_exp)*fm_W_condXZ_diffV1)
gen 	fm_deltaZ1_13 = ((((x_f_site_chi*ra_grp_exp)/sigma_v)*fm_V1_num)*(small_phi/big_phi)) +((x_f_site_chi*ra_grp_exp)*fm_W_condXZ_diffV1)
gen 	fm_deltaZ1_14 = ((((x_f_site_la*ra_grp_exp)/sigma_v)*fm_V1_num)*(small_phi/big_phi))  +((x_f_site_la*ra_grp_exp)*fm_W_condXZ_diffV1)
gen 	fm_deltaZ1_15 = ((((x_f_site_ny*ra_grp_exp)/sigma_v)*fm_V1_num)*(small_phi/big_phi))  +((x_f_site_ny*ra_grp_exp)*fm_W_condXZ_diffV1)

gen 	fm_deltaZ1_21 = ((((x_f_site_balt*ra_grp_s8)/sigma_v)*fm_V1_num)*(small_phi/big_phi)) +((x_f_site_balt*ra_grp_s8)*fm_W_condXZ_diffV1)
gen 	fm_deltaZ1_22 = ((((x_f_site_bos*ra_grp_s8)/sigma_v)*fm_V1_num)*(small_phi/big_phi))  +((x_f_site_bos*ra_grp_s8)*fm_W_condXZ_diffV1)
gen 	fm_deltaZ1_23 = ((((x_f_site_chi*ra_grp_s8)/sigma_v)*fm_V1_num)*(small_phi/big_phi))  +((x_f_site_chi*ra_grp_s8)*fm_W_condXZ_diffV1)
gen 	fm_deltaZ1_24 = ((((x_f_site_la*ra_grp_s8)/sigma_v)*fm_V1_num)*(small_phi/big_phi))   +((x_f_site_la*ra_grp_s8)*fm_W_condXZ_diffV1)
gen 	fm_deltaZ1_25 = ((((x_f_site_ny*ra_grp_s8)/sigma_v)*fm_V1_num)*(small_phi/big_phi))   +((x_f_site_ny*ra_grp_s8)*fm_W_condXZ_diffV1)

gen 	fm_deltaZ1_31 = ((((x_f_site_balt*ra_grp_control)/sigma_v)*fm_V1_num)*(small_phi/big_phi))+((x_f_site_balt*ra_grp_control)*fm_W_condXZ_diffV1)
gen 	fm_deltaZ1_32 = ((((x_f_site_bos*ra_grp_control)/sigma_v)*fm_V1_num)*(small_phi/big_phi)) +((x_f_site_bos*ra_grp_control)*fm_W_condXZ_diffV1)
gen 	fm_deltaZ1_33 = ((((x_f_site_chi*ra_grp_control)/sigma_v)*fm_V1_num)*(small_phi/big_phi)) +((x_f_site_chi*ra_grp_control)*fm_W_condXZ_diffV1)
gen 	fm_deltaZ1_34 = ((((x_f_site_la*ra_grp_control)/sigma_v)*fm_V1_num)*(small_phi/big_phi))  +((x_f_site_la*ra_grp_control)*fm_W_condXZ_diffV1)

scalar 	fm_V2_num			= ((([rho_uv2]_cons*exp([lnvar_v1]_cons))-([rho_uv1]_cons*[cov_v1v2]_cons))/det)
gen 	fm_W_condXZ_diffV2	= ((V2*exp([lnvar_v1]_cons))-(V1*[cov_v1v2]_cons))/det

gen 	fm_alpha2 	= (((1/sigma_v)*(fm_V2_num))*(small_phi/big_phi))+fm_W_condXZ_diffV2

gen 	fm_deltaZ2_11 = ((((x_f_site_balt*ra_grp_exp)/sigma_v)*fm_V2_num)*(small_phi/big_phi))+((x_f_site_balt*ra_grp_exp)*fm_W_condXZ_diffV2)
gen 	fm_deltaZ2_12 = ((((x_f_site_bos*ra_grp_exp)/sigma_v)*fm_V2_num)*(small_phi/big_phi)) +((x_f_site_bos*ra_grp_exp)*fm_W_condXZ_diffV2)
gen 	fm_deltaZ2_13 = ((((x_f_site_chi*ra_grp_exp)/sigma_v)*fm_V2_num)*(small_phi/big_phi)) +((x_f_site_chi*ra_grp_exp)*fm_W_condXZ_diffV2)
gen 	fm_deltaZ2_14 = ((((x_f_site_la*ra_grp_exp)/sigma_v)*fm_V2_num)*(small_phi/big_phi))  +((x_f_site_la*ra_grp_exp)*fm_W_condXZ_diffV2)
gen 	fm_deltaZ2_15 = ((((x_f_site_ny*ra_grp_exp)/sigma_v)*fm_V2_num)*(small_phi/big_phi))  +((x_f_site_ny*ra_grp_exp)*fm_W_condXZ_diffV2)

gen 	fm_deltaZ2_21 = ((((x_f_site_balt*ra_grp_s8)/sigma_v)*fm_V2_num)*(small_phi/big_phi)) +((x_f_site_balt*ra_grp_s8)*fm_W_condXZ_diffV2)
gen 	fm_deltaZ2_22 = ((((x_f_site_bos*ra_grp_s8)/sigma_v)*fm_V2_num)*(small_phi/big_phi))  +((x_f_site_bos*ra_grp_s8)*fm_W_condXZ_diffV2)
gen 	fm_deltaZ2_23 = ((((x_f_site_chi*ra_grp_s8)/sigma_v)*fm_V2_num)*(small_phi/big_phi))  +((x_f_site_chi*ra_grp_s8)*fm_W_condXZ_diffV2)
gen 	fm_deltaZ2_24 = ((((x_f_site_la*ra_grp_s8)/sigma_v)*fm_V2_num)*(small_phi/big_phi))   +((x_f_site_la*ra_grp_s8)*fm_W_condXZ_diffV2)
gen 	fm_deltaZ2_25 = ((((x_f_site_ny*ra_grp_s8)/sigma_v)*fm_V2_num)*(small_phi/big_phi))   +((x_f_site_ny*ra_grp_s8)*fm_W_condXZ_diffV2)

gen 	fm_deltaZ2_31 = ((((x_f_site_balt*ra_grp_control)/sigma_v)*fm_V2_num)*(small_phi/big_phi))+((x_f_site_balt*ra_grp_control)*fm_W_condXZ_diffV2)
gen 	fm_deltaZ2_32 = ((((x_f_site_bos*ra_grp_control)/sigma_v)*fm_V2_num)*(small_phi/big_phi)) +((x_f_site_bos*ra_grp_control)*fm_W_condXZ_diffV2)
gen 	fm_deltaZ2_33 = ((((x_f_site_chi*ra_grp_control)/sigma_v)*fm_V2_num)*(small_phi/big_phi)) +((x_f_site_chi*ra_grp_control)*fm_W_condXZ_diffV2)
gen 	fm_deltaZ2_34 = ((((x_f_site_la*ra_grp_control)/sigma_v)*fm_V2_num)*(small_phi/big_phi))  +((x_f_site_la*ra_grp_control)*fm_W_condXZ_diffV2)

gen 	fm_rho1_num1  = (-(V1*exp([lnvar_v2]_cons))+(V2*[cov_v1v2]_cons))/det
scalar 	fm_rho1_num2  = ((-[rho_uv1]_cons*exp([lnvar_v2]_cons))+([rho_uv2]_cons*[cov_v1v2]_cons))/det

gen X1_op0 = ((sigma_v*fm_rho1_num1)-((op0/sigma_v)*fm_rho1_num2))/(sigma_v^2)
gen X1_op1 = ((sigma_v*fm_rho1_num1)-((op1/sigma_v)*fm_rho1_num2))/(sigma_v^2)
gen X1_op2 = ((sigma_v*fm_rho1_num1)-((op2/sigma_v)*fm_rho1_num2))/(sigma_v^2)
gen X1_op3 = ((sigma_v*fm_rho1_num1)-((op3/sigma_v)*fm_rho1_num2))/(sigma_v^2)

gen 	fm_rho1 = ((normalden(op1/sigma_v)*X1_op1)-((normalden(op0/sigma_v)*X1_op0)))/big_phi if happy_scale012 == 0	
replace fm_rho1 = ((normalden(op2/sigma_v)*X1_op2)-((normalden(op1/sigma_v)*X1_op1)))/big_phi if happy_scale012 == 1	
replace fm_rho1 = ((normalden(op3/sigma_v)*X1_op3)-((normalden(op2/sigma_v)*X1_op2)))/big_phi if happy_scale012 == 2	

gen 	fm_rho2_num1 = (-(V2*exp([lnvar_v1]_cons))+(V1*[cov_v1v2]_cons))/det
scalar 	fm_rho2_num2 = (-([rho_uv2]_cons*exp([lnvar_v1]_cons))+([rho_uv1]_cons*[cov_v1v2]_cons))/det

gen X2_op0 = ((sigma_v*fm_rho2_num1)-((op0/sigma_v)*fm_rho2_num2))/(sigma_v^2)
gen X2_op1 = ((sigma_v*fm_rho2_num1)-((op1/sigma_v)*fm_rho2_num2))/(sigma_v^2)
gen X2_op2 = ((sigma_v*fm_rho2_num1)-((op2/sigma_v)*fm_rho2_num2))/(sigma_v^2)
gen X2_op3 = ((sigma_v*fm_rho2_num1)-((op3/sigma_v)*fm_rho2_num2))/(sigma_v^2)

gen 	fm_rho2 = ((normalden(op1/sigma_v)*X2_op1)-((normalden(op0/sigma_v)*X2_op0)))/big_phi if happy_scale012 == 0	
replace fm_rho2 = ((normalden(op2/sigma_v)*X2_op2)-((normalden(op1/sigma_v)*X2_op1)))/big_phi if happy_scale012 == 1	
replace fm_rho2 = ((normalden(op3/sigma_v)*X2_op3)-((normalden(op2/sigma_v)*X2_op2)))/big_phi if happy_scale012 == 2	

gen Y13	= -(((exp([lnvar_v1]_cons)^0.5)*exp([lnvar_v2]_cons))/det)-(((((((2*V1)*V2)*(exp([lnvar_v1]_cons)^0.5))*exp([lnvar_v2]_cons))*[cov_v1v2]_cons)-(((V1^2)*(exp([lnvar_v1]_cons)^0.5))*(exp([lnvar_v2]_cons)^2))-(((V2^2)*(exp([lnvar_v1]_cons)^0.5))*([cov_v1v2]_cons^2)))/(det^2))

gen 	fm_sigma1_c11  = (((((-2*V2)*[rho_uv2]_cons)*(exp([lnvar_v1]_cons)^0.5))*det)+(((2*(exp([lnvar_v1]_cons)^0.5))*exp([lnvar_v2]_cons))*mu_v_num))/(det^2)
scalar 	fm_sigma1_c12  = (((([rho_uv2]_cons^2)*(exp([lnvar_v1]_cons)^0.5))*([cov_v1v2]_cons^2))+((([rho_uv1]_cons^2)*(exp([lnvar_v1]_cons)^0.5))*(exp([lnvar_v2]_cons)^2))-/*
							*/(((((2*[rho_uv1]_cons)*[rho_uv2]_cons)*(exp([lnvar_v1]_cons)^0.5))*exp([lnvar_v2]_cons))*[cov_v1v2]_cons))/(det^2)

gen Y1_op0 = ((sigma_v*fm_sigma1_c11)-((op0/sigma_v)*fm_sigma1_c12))/(sigma_v^2)
gen Y1_op1 = ((sigma_v*fm_sigma1_c11)-((op1/sigma_v)*fm_sigma1_c12))/(sigma_v^2)
gen Y1_op2 = ((sigma_v*fm_sigma1_c11)-((op2/sigma_v)*fm_sigma1_c12))/(sigma_v^2)
gen Y1_op3 = ((sigma_v*fm_sigma1_c11)-((op3/sigma_v)*fm_sigma1_c12))/(sigma_v^2)

gen 	fm_sigma1 	= (((normalden(op1/sigma_v)*Y1_op1)-((normalden(op0/sigma_v)*Y1_op0)))/big_phi)+Y13 if happy_scale012 == 0	
replace fm_sigma1 	= (((normalden(op2/sigma_v)*Y1_op2)-((normalden(op1/sigma_v)*Y1_op1)))/big_phi)+Y13 if happy_scale012 == 1	
replace fm_sigma1 	= (((normalden(op3/sigma_v)*Y1_op3)-((normalden(op2/sigma_v)*Y1_op2)))/big_phi)+Y13 if happy_scale012 == 2	

gen Y23	= -(((exp([lnvar_v2]_cons)^0.5)*exp([lnvar_v1]_cons))/det)-(((((((2*V1)*V2)*(exp([lnvar_v2]_cons)^0.5))*exp([lnvar_v1]_cons))*[cov_v1v2]_cons)-(((V2^2)*(exp([lnvar_v2]_cons)^0.5))*(exp([lnvar_v1]_cons)^2))-(((V1^2)*(exp([lnvar_v2]_cons)^0.5))*([cov_v1v2]_cons^2)))/(det^2))

gen 	fm_sigma2_c21  = (((((-2*V1)*[rho_uv1]_cons)*(exp([lnvar_v2]_cons)^0.5))*det)+(((2*(exp([lnvar_v2]_cons)^0.5))*exp([lnvar_v1]_cons))*mu_v_num))/(det^2)
scalar 	fm_sigma2_c22  = (((([rho_uv1]_cons^2)*(exp([lnvar_v2]_cons)^0.5))*([cov_v1v2]_cons^2))+((([rho_uv2]_cons^2)*(exp([lnvar_v2]_cons)^0.5))*(exp([lnvar_v1]_cons)^2))-/*
							*/(((((2*[rho_uv1]_cons)*[rho_uv2]_cons)*(exp([lnvar_v2]_cons)^0.5))*exp([lnvar_v1]_cons))*[cov_v1v2]_cons))/(det^2)

gen Y2_op0 = ((sigma_v*fm_sigma2_c21)-((op0/sigma_v)*fm_sigma2_c22))/(sigma_v^2)
gen Y2_op1 = ((sigma_v*fm_sigma2_c21)-((op1/sigma_v)*fm_sigma2_c22))/(sigma_v^2)
gen Y2_op2 = ((sigma_v*fm_sigma2_c21)-((op2/sigma_v)*fm_sigma2_c22))/(sigma_v^2)
gen Y2_op3 = ((sigma_v*fm_sigma2_c21)-((op3/sigma_v)*fm_sigma2_c22))/(sigma_v^2)

gen 	fm_sigma2 	= (((normalden(op1/sigma_v)*Y2_op1)-((normalden(op0/sigma_v)*Y2_op0)))/big_phi)+Y23 if happy_scale012 == 0	
replace fm_sigma2 	= (((normalden(op2/sigma_v)*Y2_op2)-((normalden(op1/sigma_v)*Y2_op1)))/big_phi)+Y23 if happy_scale012 == 1	
replace fm_sigma2 	= (((normalden(op3/sigma_v)*Y2_op3)-((normalden(op2/sigma_v)*Y2_op2)))/big_phi)+Y23 if happy_scale012 == 2	

gen Z3	= ([cov_v1v2]_cons/det)+(((((V1*V2)*exp([lnvar_v2]_cons))*exp([lnvar_v1]_cons))+((V1*V2)*([cov_v1v2]_cons^2))-(((V1^2)*exp([lnvar_v2]_cons))*[cov_v1v2]_cons)-(((V2^2)*exp([lnvar_v1]_cons))*[cov_v1v2]_cons))/(det^2))

gen 	fm_sigma12_c31 = ((((V1*[rho_uv2]_cons)+(V2*[rho_uv1]_cons))*det)-((2*[cov_v1v2]_cons)*mu_v_num))/(det^2)
scalar 	fm_sigma12_c32 = (((([rho_uv1]_cons*[rho_uv2]_cons)*(exp([lnvar_v1]_cons)))*exp([lnvar_v2]_cons))+(([rho_uv1]_cons*[rho_uv2]_cons)*([cov_v1v2]_cons^2))-/*
							*/((([rho_uv1]_cons^2)*(exp([lnvar_v2]_cons)))*([cov_v1v2]_cons))-((((([rho_uv2]_cons^2)*(exp([lnvar_v1]_cons)))*[cov_v1v2]_cons))))/(det^2)

gen Z_op0 = ((sigma_v*fm_sigma12_c31)-((op0/sigma_v)*fm_sigma12_c32))/(sigma_v^2)
gen Z_op1 = ((sigma_v*fm_sigma12_c31)-((op1/sigma_v)*fm_sigma12_c32))/(sigma_v^2)
gen Z_op2 = ((sigma_v*fm_sigma12_c31)-((op2/sigma_v)*fm_sigma12_c32))/(sigma_v^2)
gen Z_op3 = ((sigma_v*fm_sigma12_c31)-((op3/sigma_v)*fm_sigma12_c32))/(sigma_v^2)

gen 	fm_sigma12 	= (((normalden(op1/sigma_v)*Z_op1)-((normalden(op0/sigma_v)*Z_op0)))/big_phi)+Z3 if happy_scale012 == 0	
replace fm_sigma12 	= (((normalden(op2/sigma_v)*Z_op2)-((normalden(op1/sigma_v)*Z_op1)))/big_phi)+Z3 if happy_scale012 == 1	
replace fm_sigma12 	= (((normalden(op3/sigma_v)*Z_op3)-((normalden(op2/sigma_v)*Z_op2)))/big_phi)+Z3 if happy_scale012 == 2	

gen N1 = (-f_c9010t_perpov_dw_z/sigma_v)*small_phi 
gen N2 = (-f_c9010t_pminorty_dw_z/sigma_v)*small_phi 

gen diff_N1_inner 		= ((normalden(op1/sigma_v))*(op1/sigma_v)) - ((normalden(op0/sigma_v))*(op0/sigma_v)) if happy_scale012 == 0									
replace diff_N1_inner 	= ((normalden(op2/sigma_v))*(op2/sigma_v)) - ((normalden(op1/sigma_v))*(op1/sigma_v)) if happy_scale012 == 1
replace diff_N1_inner 	= ((normalden(op3/sigma_v))*(op3/sigma_v)) - ((normalden(op2/sigma_v))*(op2/sigma_v)) if happy_scale012 == 2

gen diff_N1_beta = (-(f_c9010t_perpov_dw_z^2)/(sigma_v^2))*(diff_N1_inner) 

gen fm_beta1_beta1 = ((diff_N1_beta*big_phi)-(N1^2))/(big_phi^2)

gen diff_N2_inner 		= ((normalden(op1/sigma_v))*(op1/sigma_v)) - ((normalden(op0/sigma_v))*(op0/sigma_v)) if happy_scale012 == 0									
replace diff_N2_inner 	= ((normalden(op2/sigma_v))*(op2/sigma_v)) - ((normalden(op1/sigma_v))*(op1/sigma_v)) if happy_scale012 == 1
replace diff_N2_inner 	= ((normalden(op3/sigma_v))*(op3/sigma_v)) - ((normalden(op2/sigma_v))*(op2/sigma_v)) if happy_scale012 == 2

gen diff_N2_beta = (-(f_c9010t_pminorty_dw_z^2)/(sigma_v^2))*(diff_N2_inner) 

gen fm_beta2_beta2 = ((diff_N2_beta*big_phi)-(N2^2))/(big_phi^2)

replace A_balt = -(x_f_site_balt/(sigma_v))*(small_phi)
replace diff_A_balt = (-(x_f_site_balt^2)/(sigma_v^2))*(diff_N1_inner) 
gen fm_balt_balt = ((diff_A_balt*big_phi)-(A_balt^2))/(big_phi^2)

replace A_bos = -(x_f_site_bos/(sigma_v))*(small_phi)
replace diff_A_bos = (-(x_f_site_bos^2)/(sigma_v^2))*(diff_N1_inner) 
gen fm_bos_bos = ((diff_A_bos*big_phi)-(A_bos^2))/(big_phi^2)

replace A_chi = -(x_f_site_chi/(sigma_v))*(small_phi)
replace diff_A_chi = (-(x_f_site_chi^2)/(sigma_v^2))*(diff_N1_inner) 
gen fm_chi_chi = ((diff_A_chi*big_phi)-(A_chi^2))/(big_phi^2)

replace A_la = -(x_f_site_la/(sigma_v))*(small_phi)
replace diff_A_la = (-(x_f_site_la^2)/(sigma_v^2))*(diff_N1_inner) 
gen fm_la_la = ((diff_A_la*big_phi)-(A_la^2))/(big_phi^2)

gen diff_N1_balt 	= (-(f_c9010t_perpov_dw_z*x_f_site_balt)/(sigma_v^2))*(diff_N1_inner) 
gen diff_N1_bos 	= (-(f_c9010t_perpov_dw_z*x_f_site_bos)/(sigma_v^2))*(diff_N1_inner) 
gen diff_N1_chi 	= (-(f_c9010t_perpov_dw_z*x_f_site_chi)/(sigma_v^2))*(diff_N1_inner) 
gen diff_N1_la 		= (-(f_c9010t_perpov_dw_z*x_f_site_la)/(sigma_v^2))*(diff_N1_inner) 

gen diff_N2_balt 	= (-(f_c9010t_pminorty_dw_z*x_f_site_balt)/(sigma_v^2))*(diff_N2_inner) 
gen diff_N2_bos 	= (-(f_c9010t_pminorty_dw_z*x_f_site_bos)/(sigma_v^2))*(diff_N2_inner) 
gen diff_N2_chi 	= (-(f_c9010t_pminorty_dw_z*x_f_site_chi)/(sigma_v^2))*(diff_N2_inner) 
gen diff_N2_la 		= (-(f_c9010t_pminorty_dw_z*x_f_site_la)/(sigma_v^2))*(diff_N2_inner) 

replace diff_D_balt = -(x_f_site_balt/(sigma_v))*(small_phi)
replace diff_D_bos 	= -(x_f_site_bos/(sigma_v))*(small_phi)
replace diff_D_chi 	= -(x_f_site_chi/(sigma_v))*(small_phi)
replace diff_D_la 	= -(x_f_site_la/(sigma_v))*(small_phi)

gen fm_balt_beta1 	= ((diff_N1_balt*big_phi)-(N1*diff_D_balt))/(big_phi^2)
gen fm_bos_beta1 	= ((diff_N1_bos*big_phi)-(N1*diff_D_bos))/(big_phi^2)
gen fm_chi_beta1 	= ((diff_N1_chi*big_phi)-(N1*diff_D_chi))/(big_phi^2)
gen fm_la_beta1 	= ((diff_N1_la*big_phi)-(N1*diff_D_la))/(big_phi^2)

gen fm_balt_beta2 	= ((diff_N2_balt*big_phi)-(N2*diff_D_balt))/(big_phi^2)
gen fm_bos_beta2 	= ((diff_N2_bos*big_phi)-(N2*diff_D_bos))/(big_phi^2)
gen fm_chi_beta2 	= ((diff_N2_chi*big_phi)-(N2*diff_D_chi))/(big_phi^2)
gen fm_la_beta2 	= ((diff_N2_la*big_phi)-(N2*diff_D_la))/(big_phi^2)

gen combine_beta1_beta1 	= (fm_beta1*fm_beta1) + fm_beta1_beta1 
gen combine_beta2_beta2 	= (fm_beta2*fm_beta2) + fm_beta2_beta2 
replace combine_balt_balt 	= (fm_balt*fm_balt) + fm_balt_balt 
replace combine_bos_bos 	= (fm_bos*fm_bos) + fm_bos_bos 
replace combine_chi_chi 	= (fm_chi*fm_chi) + fm_chi_chi 
replace combine_la_la 		= (fm_la*fm_la) + fm_la_la 

gen combine_beta1_balt 	= (fm_beta1*fm_balt) + fm_balt_beta1 
gen combine_beta1_bos 	= (fm_beta1*fm_bos) + fm_bos_beta1 
gen combine_beta1_chi 	= (fm_beta1*fm_chi) + fm_chi_beta1 
gen combine_beta1_la 	= (fm_beta1*fm_la) + fm_la_beta1 

gen combine_beta2_balt 	= (fm_beta2*fm_balt) + fm_balt_beta2 
gen combine_beta2_bos 	= (fm_beta2*fm_bos) + fm_bos_beta2 
gen combine_beta2_chi 	= (fm_beta2*fm_chi) + fm_chi_beta2 
gen combine_beta2_la 	= (fm_beta2*fm_la) + fm_la_beta2 

* now run regressions of 1st and 2nd derivatives on column of ones to get the IM test statistic  

regress ones fm_beta1 fm_beta2 fm_balt fm_bos fm_chi fm_la fm_c2 fm_c1 fm_alpha1 fm_deltaZ1_11 fm_deltaZ1_12 fm_deltaZ1_13 fm_deltaZ1_14 fm_deltaZ1_15 /*
			*/fm_deltaZ1_21 fm_deltaZ1_22 fm_deltaZ1_23 fm_deltaZ1_24 fm_deltaZ1_25 fm_deltaZ1_31 fm_deltaZ1_32 fm_deltaZ1_33 fm_deltaZ1_34 fm_alpha2 /*
			*/fm_deltaZ2_11 fm_deltaZ2_12 fm_deltaZ2_13 fm_deltaZ2_14 fm_deltaZ2_15 fm_deltaZ2_21 fm_deltaZ2_22 fm_deltaZ2_23 fm_deltaZ2_24 fm_deltaZ2_25 /*
			*/fm_deltaZ2_31 fm_deltaZ2_32 fm_deltaZ2_33 fm_deltaZ2_34 fm_rho1 fm_rho2 fm_sigma1 fm_sigma2 fm_sigma12, noconstant

* n=3263 and R-square=0.0036 in the regression below which gives an IM test statistic of 11.7468. The p-value reported in the paper comes from 
* a chi-square distribution with 3 dof and equals 0.0083  

regress ones fm_beta1 fm_beta2 fm_balt fm_bos fm_chi fm_la fm_c2 fm_c1 fm_alpha1 fm_deltaZ1_11 fm_deltaZ1_12 fm_deltaZ1_13 fm_deltaZ1_14 fm_deltaZ1_15 /*
			*/fm_deltaZ1_21 fm_deltaZ1_22 fm_deltaZ1_23 fm_deltaZ1_24 fm_deltaZ1_25 fm_deltaZ1_31 fm_deltaZ1_32 fm_deltaZ1_33 fm_deltaZ1_34 fm_alpha2 /*
			*/fm_deltaZ2_11 fm_deltaZ2_12 fm_deltaZ2_13 fm_deltaZ2_14 fm_deltaZ2_15 fm_deltaZ2_21 fm_deltaZ2_22 fm_deltaZ2_23 fm_deltaZ2_24 fm_deltaZ2_25 /*
			*/fm_deltaZ2_31 fm_deltaZ2_32 fm_deltaZ2_33 fm_deltaZ2_34 fm_rho1 fm_rho2 fm_sigma1 fm_sigma2 fm_sigma12 combine_beta1_beta1 combine_balt_balt combine_beta1_balt, noconstant
