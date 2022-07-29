
** LINEAR MODELS, APPENDIX C

** TABLE C1 

** Panel A, linear models assuming exogeneity

regress happy_scale012 $site_covs f_c9010t_perpov_dw_z [pw=$wt]
est store ne_l_pov_1

regress happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z [pw=$wt]
est store ne_l_pov_2

regress happy_scale012 $site_covs f_c9010t_pminorty_dw_z [pw=$wt]
est store ne_l_min_1

regress happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_pminorty_dw_z [pw=$wt]
est store ne_l_min_2

regress happy_scale012 $site_covs f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z [pw=$wt]
est store ne_l_povmin_1

regress happy_scale012 $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z [pw=$wt]
est store ne_l_povmin_2

estout ne_l_pov_1 ne_l_min_1 ne_l_povmin_1 ne_l_pov_2 ne_l_min_2 ne_l_povmin_2 ///
using "${pathres}tableC1_${date}.tex", cells("b(star label() fmt(4))" "se(par label() fmt(4))") ///
replace mlabel(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) stardetach label ///
varwidth(12) modelwidth(10) style(tex) keep(f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z) stats(N, fmt(%9.0f)) ///
varlabels(f_c9010t_perpov_dw_z "$ \beta_{Poverty} $" f_c9010t_pminorty_dw_z "$ \beta_{Minority} $") ///
title(Linear model estimation (OLS and IV) of neighborhood effects on SWB \label{ne_l}) ///
prehead("\begin{table}[!htbp] \begin{center} \begin{threeparttable} \caption{@title}\renewcommand{\arraystretch}{1.2} \small " ///
"\begin{tabular}{@{}l r@{}l r@{}l r@{}l r@{}l r@{}l r@{}l} \toprule " ///
" & \multicolumn{2}{c}{(1)} & \multicolumn{2}{c}{(2)} & \multicolumn{2}{c}{(3)} & \multicolumn{2}{c}{(4)} & \multicolumn{2}{c}{(5)} & \multicolumn{2}{c}{(6)}\\" ///
" \midrule ") ///
posthead("\textbf{Panel A: OLS estimation} & & & & & & & & & & & \\") postfoot("\midrule") 

** Panel B, linear models but using IV

ivregress 2sls happy_scale012 $site_covs (f_c9010t_perpov_dw_z = i.ra_group#i.ra_site) [pw=$wt]
est store ne_liv_pov_1

ivregress 2sls happy_scale012 $site_covs x_f_release1 $x_covariates (f_c9010t_perpov_dw_z = i.ra_group#i.ra_site) [pw=$wt]
est store ne_liv_pov_2

ivregress 2sls happy_scale012 $site_covs (f_c9010t_pminorty_dw_z = i.ra_group#i.ra_site) [pw=$wt]
est store ne_liv_min_1

ivregress 2sls happy_scale012 $site_covs x_f_release1 $x_covariates (f_c9010t_pminorty_dw_z = i.ra_group#i.ra_site) [pw=$wt]
est store ne_liv_min_2

ivregress 2sls happy_scale012 $site_covs (f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z = i.ra_group#i.ra_site) [pw=$wt]
est store ne_liv_povmin_1

ivregress 2sls happy_scale012 $site_covs x_f_release1 $x_covariates (f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z = i.ra_group#i.ra_site) [pw=$wt]
est store ne_liv_povmin_2

estout ne_liv_pov_1 ne_liv_min_1 ne_liv_povmin_1 ne_liv_pov_2 ne_liv_min_2 ne_liv_povmin_2 ///
using "${pathres}tableC1_${date}.tex", cells("b(star label() fmt(4))" "se(par label() fmt(4))") ///
append mlabel(none) collabels(none) ///
starlevels(* 0.10 ** 0.05 *** 0.01) stardetach label ///
varwidth(12) modelwidth(10) style(tex) keep(f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z) stats(N, fmt(%9.0f)) ///
varlabels(f_c9010t_perpov_dw_z "$ \beta_{Poverty} $" f_c9010t_pminorty_dw_z "$ \beta_{Minority} $") ///
posthead("\textbf{Panel B: IV estimation} & & & & & & & & & & & \\") prefoot(" ") ///
postfoot("\bottomrule \end{tabular}" ///
"\begin{TableNotes} \scriptsize" ///
"\item \textit{Notes:} The dependent variable is Subjective Well Being (SWB) which takes the value zero for not too happy, one for pretty happy and two for very happy; columns (1)-(3) use a set of dummy variables for randomization site as covariates X while columns (4)-(6) use a complete set of baseline characteristics (as given in Table \ref{tab:table1}), and whether a sample adult was included in the first release of the long-term evaluation survey fielding period, as covariates X; all regressions are weighted; * p-value $<0.10$, ** p-value $<0.05$, *** p-value $<0.01$." ///
"\item \textit{Source:} Data from ICPSR Study 34860: Moving to Opportunity: Final Impacts Evaluation Science Article Data, 2008-2010." ///
"\end{TableNotes} \end{threeparttable} \end{center} \end{table}")

** TABLE C2

gen ra_exp_exps8 = ra_poolgrp_exps8 if ra_group==1|ra_group==3
gen ra_sec8_exps8 = ra_poolgrp_exps8 if ra_group==2|ra_group==3

regress happy_scale012 ra_poolgrp_exps8 $site_covs [pw=$wt]
est store itt_l_pool_1

regress happy_scale012 ra_poolgrp_exps8 $site_covs x_f_release1 $x_covariates [pw=$wt]
est store itt_l_pool_2

regress happy_scale012 ra_exp_exps8 $site_covs [pw=$wt]
est store itt_l_exp_1

regress happy_scale012 ra_exp_exps8 $site_covs x_f_release1 $x_covariates [pw=$wt]
est store itt_l_exp_2

regress happy_scale012 ra_sec8_exps8 $site_covs [pw=$wt]
est store itt_l_sec8_1

regress happy_scale012 ra_sec8_exps8 $site_covs x_f_release1 $x_covariates [pw=$wt]
est store itt_l_sec8_2

estout itt_l_pool_1 itt_l_exp_1 itt_l_sec8_1 itt_l_pool_2 itt_l_exp_2 itt_l_sec8_2 ///
using "${pathres}tableC2_${date}.tex", cells("b(star label() fmt(4))" "se(par label() fmt(4))") ///
replace mlabel(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) stardetach label ///
varwidth(12) modelwidth(10) style(tex) keep(ra_poolgrp_exps8 ra_exp_exps8 ra_sec8_exps8) stats(N, fmt(%9.0f)) ///
varlabels(ra_poolgrp_exps8 "$ Z= $ Any MTO voucher" ra_exp_exps8 "$ Z= $ MTO low poverty voucher" ra_sec8_exps8 "$ Z= $ MTO section 8 voucher") ///
title(Linear model estimation (ITT) of neighborhood effects on SWB \label{itt_l}) ///
prehead("\begin{table}[!htbp] \begin{center} \begin{threeparttable} \caption{@title}\renewcommand{\arraystretch}{1.2} \small " ///
"\begin{tabular}{@{}l r@{}l r@{}l r@{}l r@{}l r@{}l r@{}l} \toprule " ///
" & \multicolumn{2}{c}{(1)} & \multicolumn{2}{c}{(2)} & \multicolumn{2}{c}{(3)} & \multicolumn{2}{c}{(4)} & \multicolumn{2}{c}{(5)} & \multicolumn{2}{c}{(6)}\\" ///
" \midrule ") ///
posthead(" ") prefoot("\midrule ") ///
postfoot("\bottomrule \end{tabular}" ///
"\begin{TableNotes} \scriptsize" ///
"\item \textit{Notes:} The dependent variable is Subjective Well Being (SWB) which takes the value zero for not too happy, one for pretty happy and two for very happy; columns (1)-(3) use a set of dummy variables for randomization site as covariates X while columns (4)-(6) use a complete set of baseline characteristics (as given in Table \ref{tab:table1}), and whether a sample adult was included in the first release of the long-term evaluation survey fielding period, as covariates X; all regressions are weighted; * p-value $<0.10$, ** p-value $<0.05$, *** p-value $<0.01$." ///
"\item \textit{Source:} Data from ICPSR Study 34860: Moving to Opportunity: Final Impacts Evaluation Science Article Data, 2008-2010." ///
"\end{TableNotes} \end{threeparttable} \end{center} \end{table}")
