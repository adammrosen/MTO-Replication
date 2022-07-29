
** TABLE 1, but ordered logit instead of ordered probit 

ologit happy_scale012 $site_covs f_c9010t_perpov_dw_z [pw=$wt]
est store ol_ne_nl_pov_1

ologit happy_scale012 $site_covs f_c9010t_pminorty_dw_z [pw=$wt]
est store ol_ne_nl_min_1

ologit happy_scale012 $site_covs f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z [pw=$wt]
est store ol_ne_nl_povmin_1

ologit happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z [pw=$wt]
est store ol_ne_nl_pov_2

ologit happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_pminorty_dw_z [pw=$wt]
est store ol_ne_nl_min_2

ologit happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z [pw=$wt]
est store ol_ne_nl_povmin_2

estout ol_ne_nl_pov_1 ol_ne_nl_min_1 ol_ne_nl_povmin_1 ol_ne_nl_pov_2 ol_ne_nl_min_2 ol_ne_nl_povmin_2 ///
using "${pathres}_${date}_ol_ne_nl_w.tex", cells("b(star label() fmt(4))" "se(par label() fmt(4))") ///
replace mlabel(none) collabels(none) ///
starlevels(* 0.10 ** 0.05 *** 0.01) stardetach label ///
varwidth(12) modelwidth(10) style(tex) ///
drop(tmp_female x_rad_ad_le_35 x_rad_ad_36_40 x_rad_ad_41_45 x_rad_ad_46_50 ///
x_rad_ad_ethrace_black_nh x_rad_ad_ethrace_hisp ///
x_f_ad_nevmarr x_f_ad_parentu18 x_f_ad_working x_f_ad_edinsch x_f_ad_edgradhs x_f_ad_edged x_f_hh_afdc ///
rad_svy_bl_totincm_2009d ///
x_f_hh_car x_f_hh_disabl x_f_hh_noteens x_f_hh_size2 x_f_hh_size3 x_f_hh_size4 ///
x_f_hh_victim x_f_hood_unsafenit x_f_hood_verydissat x_f_hood_5y x_f_hous_mov3tm  ///
x_f_hood_nofamily x_f_hood_nofriend x_f_hood_chat x_f_hood_nbrkid x_f_hous_fndapt x_f_hous_sec8bef ///
x_f_hous_movdrgs x_f_hous_movschl cov_hous_movapt cov_hous_movjob x_f_release1) ///
order(f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z x_f_site_balt x_f_site_bos x_f_site_chi x_f_site_la) ///
rename(_cons cons) ///
varlabels(f_c9010t_perpov_dw_z "$ \beta_{Poverty} $" f_c9010t_pminorty_dw_z "$ \beta_{Minority} $" x_f_site_balt "$ \gamma_{Baltomore} $" x_f_site_bos "$ \gamma_{Boston} $" ///
x_f_site_chi "$ \gamma_{Chicago} $" x_f_site_la "$ \gamma_{LA} $" /:cut1 "$ c_{1} $" /:cut2 "$ c_{2} $") ///
stats(N, fmt(%9.0f)) ///
title(Neighborhood effects on SWB, Ordered Logit model \label{ol_ne_nl_w}) ///
prehead("\begin{table}[!htbp] \begin{center} \begin{threeparttable} \caption{@title}\renewcommand{\arraystretch}{1.2} \small " ///
"\begin{tabular}{@{}l r@{}l r@{}l r@{}l r@{}l r@{}l r@{}l} \toprule " ///
" & \multicolumn{2}{c}{(1)} & \multicolumn{2}{c}{(2)} & \multicolumn{2}{c}{(3)} & \multicolumn{2}{c}{(4)} & \multicolumn{2}{c}{(5)} & \multicolumn{2}{c}{(6)}\\" ///
" \midrule ") ///
prefoot("\midrule ") ///
postfoot("\bottomrule \end{tabular}" ///
"\begin{TableNotes} \scriptsize" ///
"\item \textit{Notes:} The dependent variable is Subjective Well Being (SWB) which takes the value zero for not too happy, one for pretty happy and two for very happy; columns (1)-(3) use a set of dummy variables for randomization site as covariates X while columns (4)-(6) use a complete set of baseline characteristics, and whether a sample adult was included in the first release of the long-term evaluation survey fielding period, as covariates X; all regressions are weighted; * p-value $<0.10$, ** p-value $<0.05$, *** p-value $<0.01$." ///
"\item \textit{Source:} Data from ICPSR Study 34860: Moving to Opportunity: Final Impacts Evaluation Science Article Data, 2008-2010." ///
"\end{TableNotes} \end{threeparttable} \end{center} \end{table}")

** ordered logit probabilities and marginal effects

** W = hood poverty, with covariates

ologit happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z [pw=$wt]

*** conditional probabilities that Y = 0,1,2

forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=(-6(1)4)) predict(outcome(`y'))
		
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(0(0.2)1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black))
	graph export ${pathres}prob_ol_y`y'_x_pov_med_ny_${date}.eps, replace
	}

*** conditional marginal effects for Y = 0,1,2

forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=(-6(1)4)) dydx(f_c9010t_perpov_dw_z) predict(outcome(`y')) 
	
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(-0.1(0.05)0.1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black)) 
	graph export ${pathres}me_ol_y`y'_x_pov_med_ny_${date}.eps, replace
	}
	
** W = hood minority, with covariates

ologit happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_pminorty_dw_z [pw=$wt]

*** conditional probabilities that Y = 0,1,2

forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=(-6(1)4)) predict(outcome(`y'))  
	
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(0(0.2)1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black))
	graph export ${pathres}prob_ol_y`y'_x_min_med_ny_${date}.eps, replace
	}

*** conditional marginal effects for Y = 0,1,2

forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=(-6(1)4)) dydx(f_c9010t_pminorty_dw_z) predict(outcome(`y'))  
	
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(-0.1(0.05)0.1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black)) 
	graph export ${pathres}me_ol_y`y'_x_min_med_ny_${date}.eps, replace
	}
	
** W = (hood poverty, hood minority), with covariates
	
ologit happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z [pw=$wt]

*** conditional probabilities that Y = 0,1,2

forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=.4965963 f_c9010t_perpov_dw_z=(-6(1)4)) predict(outcome(`y')) 
	
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(0(0.2)1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black))
	graph export ${pathres}prob_ol_y`y'_x_povmin1_med_ny_${date}.eps, replace
	}
	
forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=-.182967 f_c9010t_pminorty_dw_z=(-6(1)4)) predict(outcome(`y')) 
	
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(0(0.2)1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black))
	graph export ${pathres}prob_ol_y`y'_x_povmin2_med_ny_${date}.eps, replace
	}

*** conditional marginal effects for Y = 0,1,2

forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=.4965963 f_c9010t_perpov_dw_z=(-6(1)4)) dydx(f_c9010t_perpov_dw_z) predict(outcome(`y')) 
	
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(-0.1(0.05)0.1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black)) 
	graph export ${pathres}me_ol_y`y'_x_povmin1_med_ny_${date}.eps, replace
	}
	
forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=-.182967 f_c9010t_pminorty_dw_z=(-6(1)4)) dydx(f_c9010t_pminorty_dw_z) predict(outcome(`y')) 
	
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(-0.1(0.05)0.1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black))  
	graph export ${pathres}me_ol_y`y'_x_povmin2_med_ny_${date}.eps, replace
	}
	
estimates clear 

** TABLE 1, but multinomial logit instead of ordered probit 

mlogit happy_scale012 $site_covs f_c9010t_perpov_dw_z [pw=$wt]
est store ml_ne_nl_pov_1

mlogit happy_scale012 $site_covs f_c9010t_pminorty_dw_z [pw=$wt]
est store ml_ne_nl_min_1

mlogit happy_scale012 $site_covs f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z [pw=$wt]
est store ml_ne_nl_povmin_1

mlogit happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z [pw=$wt]
est store ml_ne_nl_pov_2

mlogit happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_pminorty_dw_z [pw=$wt]
est store ml_ne_nl_min_2

mlogit happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z [pw=$wt]
est store ml_ne_nl_povmin_2

estout ml_ne_nl_pov_1 ml_ne_nl_min_1 ml_ne_nl_povmin_1 ml_ne_nl_pov_2 ml_ne_nl_min_2 ml_ne_nl_povmin_2 ///
using "${pathres}_${date}_ml_ne_nl_w.tex", cells("b(star label() fmt(4))" "se(par label() fmt(4))") ///
replace mlabel(none) collabels(none) noomitted ///
starlevels(* 0.10 ** 0.05 *** 0.01) stardetach label ///
varwidth(12) modelwidth(10) style(tex) ///
drop(tmp_female x_rad_ad_le_35 x_rad_ad_36_40 x_rad_ad_41_45 x_rad_ad_46_50 ///
x_rad_ad_ethrace_black_nh x_rad_ad_ethrace_hisp ///
x_f_ad_nevmarr x_f_ad_parentu18 x_f_ad_working x_f_ad_edinsch x_f_ad_edgradhs x_f_ad_edged x_f_hh_afdc ///
rad_svy_bl_totincm_2009d ///
x_f_hh_car x_f_hh_disabl x_f_hh_noteens x_f_hh_size2 x_f_hh_size3 x_f_hh_size4 ///
x_f_hh_victim x_f_hood_unsafenit x_f_hood_verydissat x_f_hood_5y x_f_hous_mov3tm  ///
x_f_hood_nofamily x_f_hood_nofriend x_f_hood_chat x_f_hood_nbrkid x_f_hous_fndapt x_f_hous_sec8bef ///
x_f_hous_movdrgs x_f_hous_movschl cov_hous_movapt cov_hous_movjob x_f_release1) ///
order(f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z x_f_site_balt x_f_site_bos x_f_site_chi x_f_site_la) ///
varlabels(0:f_c9010t_perpov_dw_z "$ \beta_{Poverty} $" 0:f_c9010t_pminorty_dw_z "$ \beta_{Minority} $" 0:x_f_site_balt "$ \gamma_{Baltomore} $" 0:x_f_site_bos "$ \gamma_{Boston} $" ///
0:x_f_site_chi "$ \gamma_{Chicago} $" 0:x_f_site_la "$ \gamma_{LA} $" 0:_cons "$ c_{0} $" 2:f_c9010t_perpov_dw_z "$ \beta_{Poverty} $" 2:f_c9010t_pminorty_dw_z "$ \beta_{Minority} $" ///
2:x_f_site_balt "$ \gamma_{Baltomore} $" 2:x_f_site_bos "$ \gamma_{Boston} $" 2:x_f_site_chi "$ \gamma_{Chicago} $" 2:x_f_site_la "$ \gamma_{LA} $" 2:_cons "$ c_{2} $") ///
stats(N, fmt(%9.0f)) ///
title(Neighborhood effects on SWB, Multinomial Logit model \label{ml_ne_nl_w}) ///
prehead("\begin{table}[!htbp] \begin{center} \begin{threeparttable} \caption{@title}\renewcommand{\arraystretch}{1.2} \small " ///
"\begin{tabular}{@{}l r@{}l r@{}l r@{}l r@{}l r@{}l r@{}l} \toprule " ///
" & \multicolumn{2}{c}{(1)} & \multicolumn{2}{c}{(2)} & \multicolumn{2}{c}{(3)} & \multicolumn{2}{c}{(4)} & \multicolumn{2}{c}{(5)} & \multicolumn{2}{c}{(6)}\\" ///
" \midrule ") ///
prefoot("\midrule ") ///
postfoot("\bottomrule \end{tabular}" ///
"\begin{TableNotes} \scriptsize" ///
"\item \textit{Notes:} The dependent variable is Subjective Well Being (SWB) which takes the value zero for not too happy, one for pretty happy and two for very happy; columns (1)-(3) use a set of dummy variables for randomization site as covariates X while columns (4)-(6) use a complete set of baseline characteristics, and whether a sample adult was included in the first release of the long-term evaluation survey fielding period, as covariates X; all regressions are weighted; * p-value $<0.10$, ** p-value $<0.05$, *** p-value $<0.01$." ///
"\item \textit{Source:} Data from ICPSR Study 34860: Moving to Opportunity: Final Impacts Evaluation Science Article Data, 2008-2010." ///
"\end{TableNotes} \end{threeparttable} \end{center} \end{table}")

** multinomial logit probabilities and marginal effects

** W = hood poverty, with covariates

mlogit happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z [pw=$wt]

*** conditional probabilities that Y = 0,1,2

forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=(-6(1)4)) predict(outcome(`y'))
		
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(0(0.2)1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black))
	graph export ${pathres}prob_ml_y`y'_x_pov_med_ny_${date}.eps, replace
	}

*** conditional marginal effects for Y = 0,1,2

forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=(-6(1)4)) dydx(f_c9010t_perpov_dw_z) predict(outcome(`y')) 
	
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(-0.1(0.05)0.1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black)) 
	graph export ${pathres}me_ml_y`y'_x_pov_med_ny_${date}.eps, replace
	}
	
** W = hood minority, with covariates

mlogit happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_pminorty_dw_z [pw=$wt]

*** conditional probabilities that Y = 0,1,2

forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=(-6(1)4)) predict(outcome(`y'))  
	
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(0(0.2)1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black))
	graph export ${pathres}prob_ml_y`y'_x_min_med_ny_${date}.eps, replace
	}

*** conditional marginal effects for Y = 0,1,2

forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=(-6(1)4)) dydx(f_c9010t_pminorty_dw_z) predict(outcome(`y'))  
	
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(-0.1(0.05)0.1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black)) 
	graph export ${pathres}me_ml_y`y'_x_min_med_ny_${date}.eps, replace
	}
	
** W = (hood poverty, hood minority), with covariates
	
mlogit happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z [pw=$wt]

*** conditional probabilities that Y = 0,1,2

forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=.4965963 f_c9010t_perpov_dw_z=(-6(1)4)) predict(outcome(`y')) 
	
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(0(0.2)1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black))
	graph export ${pathres}prob_ml_y`y'_x_povmin1_med_ny_${date}.eps, replace
	}
	
forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=-.182967 f_c9010t_pminorty_dw_z=(-6(1)4)) predict(outcome(`y')) 
	
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(0(0.2)1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black))
	graph export ${pathres}prob_ml_y`y'_x_povmin2_med_ny_${date}.eps, replace
	}

*** conditional marginal effects for Y = 0,1,2

forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_pminorty_dw_z=.4965963 f_c9010t_perpov_dw_z=(-6(1)4)) dydx(f_c9010t_perpov_dw_z) predict(outcome(`y')) 
	
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(-0.1(0.05)0.1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black)) 
	graph export ${pathres}me_ml_y`y'_x_povmin1_med_ny_${date}.eps, replace
	}
	
forvalues y = 0/2 {
	margins, at(${x_fix_long} f_c9010t_perpov_dw_z=-.182967 f_c9010t_pminorty_dw_z=(-6(1)4)) dydx(f_c9010t_pminorty_dw_z) predict(outcome(`y')) 
	
	marginsplot, title(" ") xtitle(" ") ytitle(" ") ylabel(-0.1(0.05)0.1) xlabel(-6(2)4) plotopts(mcolor(black) lcolor(black)) ciopts(lcolor(black))  
	graph export ${pathres}me_ml_y`y'_x_povmin2_med_ny_${date}.eps, replace
	}
	
estimates clear 