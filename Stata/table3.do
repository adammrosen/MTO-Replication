
** TABLE 3

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs f_c9010t_perpov_dw_z, nocons) (cut1:) (lndiff:) (xb2:f_c9010t_perpov_dw_z = interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml maximize
est store pov_t_nox

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs f_c9010t_pminorty_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_pminorty_dw_z = interaction*) (lnsigma:) (rho:) [pw=$wt], /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml init xb1:x_f_site_bal=0.3148 x_f_site_bos=0.1572 x_f_site_chi=0.2667 x_f_site_la=0.1045 f_c9010t_pminorty_dw_z=-0.0588
//ml init eq1:x_f_site_bal=0.3150 x_f_site_bos=0.1577 x_f_site_chi=0.2666 x_f_site_la=0.1045 f_c9010t_pminorty_dw_z=-0.0585
ml maximize, difficult
est store min_t_nox

ml model lf myoprobit3_lf (xb1:happy_scale012 = $site_covs f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z, nocons) (cut1:) (lndiff:) /*
						*/(xb2:f_c9010t_perpov_dw_z = interaction*) (lnvar_v1:) (rho_uv1:) /*
						*/(xb3:f_c9010t_pminorty_dw_z = interaction*) (lnvar_v2:) (rho_uv2:) (cov_v1v2:) [pw=$wt], /*
						*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
						*/diparm(lnvar_v1, exp label(var_v1)) /*
						*/diparm(lnvar_v2, exp label(var_v2)) 
ml init xb1:x_f_site_bal=0.3928 x_f_site_bos=0.4820 x_f_site_chi=0.2571 x_f_site_la=0.0768 f_c9010t_perpov_dw_z=-0.2857 f_c9010t_pminorty_dw_z=0.3206 cut1:_cons=-0.4024
ml maximize, difficult
est store povmin_t_nox

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_perpov_dw_z = x_f_release1 $x_covariates interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml init xb1:x_f_site_bal=0.1825 x_f_site_bos=0.0528 x_f_site_chi=0.1811 x_f_site_la=0.0674 f_c9010t_perpov_dw_z=-0.1610
ml maximize, difficult
est store pov_t_x

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_pminorty_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_pminorty_dw_z = x_f_release1 $x_covariates interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml maximize
est store min_t_x

ml model lf myoprobit3_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z, nocons) (cut1:) (lndiff:) /*
						*/(xb2:f_c9010t_perpov_dw_z = x_f_release1 $x_covariates interaction*) (lnvar_v1:) (rho_uv1:) /*
						*/(xb3:f_c9010t_pminorty_dw_z = x_f_release1 $x_covariates interaction*) (lnvar_v2:) (rho_uv2:) (cov_v1v2:) [pw=$wt], /*
						*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
						*/diparm(lnvar_v1, exp label(var_v1)) /*
						*/diparm(lnvar_v2, exp label(var_v2)) 
ml init xb1:x_f_site_bal=0.3039 x_f_site_bos=0.3386 x_f_site_chi=0.1527 x_f_site_la=0.0639 f_c9010t_perpov_dw_z=-0.2998 f_c9010t_pminorty_dw_z=0.3151 cut1:_cons=-0.8165
ml maximize, difficult
est store povmin_t_x

estout pov_t_nox min_t_nox povmin_t_nox pov_t_x min_t_x povmin_t_x ///
using "${pathres}table3_${date}.tex", cells("b(star label() fmt(4))" "se(par label() fmt(4))") ///
replace mlabel(none) collabels(none) eqlabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) stardetach label ///
varwidth(12) modelwidth(10) style(tex) drop(${x_covariates} x_f_release1 xb2:_cons xb2:interaction* lnvar_v1:_cons xb3:_cons xb3:interaction* lnvar_v2:_cons  lnsigma:_cons) stats(N, fmt(%9.0f)) ///
order(f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z x_f_site_balt x_f_site_bos x_f_site_chi x_f_site_la cut1:_cons lndiff:_cons rho:_cons rho_uv1:_cons rho_uv2:_cons cov_v1v2:_cons) ///
varlabels(f_c9010t_perpov_dw_z "$ \beta_{Poverty} $" f_c9010t_pminorty_dw_z "$ \beta_{Minority} $" x_f_site_balt "$ \gamma_{Baltomore} $" x_f_site_bos "$ \gamma_{Boston} $" ///
x_f_site_chi "$ \gamma_{Chicago} $" x_f_site_la "$ \gamma_{LA} $" cut1:_cons "$ c_{1} $" lndiff:_cons "$ \ln(c_{2}-c_{1}) $" rho:_cons "$ \rho $" rho_uv1:_cons "$ \rho_{1} $" rho_uv2:_cons "$ \rho_{2} $" cov_v1v2:_cons "$ cov(v_{1},v_{2}) $") ///
title(Triangular IV estimation of neighborhood effects on SWB \label{tab:tiv}) ///
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

** recover parameter estimates 

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs f_c9010t_perpov_dw_z, nocons) (cut1:) (lndiff:) (xb2:f_c9010t_perpov_dw_z = interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml maximize

nlcom exp([lndiff]_cons)+[cut1]_cons

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs f_c9010t_pminorty_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_pminorty_dw_z = interaction*) (lnsigma:) (rho:) [pw=$wt], /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml init xb1:x_f_site_bal=0.3148 x_f_site_bos=0.1572 x_f_site_chi=0.2667 x_f_site_la=0.1045 f_c9010t_pminorty_dw_z=-0.0588
//ml init eq1:x_f_site_bal=0.3150 x_f_site_bos=0.1577 x_f_site_chi=0.2666 x_f_site_la=0.1045 f_c9010t_pminorty_dw_z=-0.0585
ml maximize, difficult

nlcom exp([lndiff]_cons)+[cut1]_cons

ml model lf myoprobit3_lf (xb1:happy_scale012 = $site_covs f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z, nocons) (cut1:) (lndiff:) /*
						*/(xb2:f_c9010t_perpov_dw_z = interaction*) (lnvar_v1:) (rho_uv1:) /*
						*/(xb3:f_c9010t_pminorty_dw_z = interaction*) (lnvar_v2:) (rho_uv2:) (cov_v1v2:) [pw=$wt], /*
						*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
						*/diparm(lnvar_v1, exp label(var_v1)) /*
						*/diparm(lnvar_v2, exp label(var_v2)) 
ml init xb1:x_f_site_bal=0.3928 x_f_site_bos=0.4820 x_f_site_chi=0.2571 x_f_site_la=0.0768 f_c9010t_perpov_dw_z=-0.2857 f_c9010t_pminorty_dw_z=0.3206 cut1:_cons=-0.4024
ml maximize, difficult

nlcom exp([lndiff]_cons)+[cut1]_cons
nlcom (exp([lnvar_v1]_cons))^0.5
nlcom (exp([lnvar_v2]_cons))^0.5

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_perpov_dw_z = x_f_release1 $x_covariates interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml init xb1:x_f_site_bal=0.1825 x_f_site_bos=0.0528 x_f_site_chi=0.1811 x_f_site_la=0.0674 f_c9010t_perpov_dw_z=-0.1610
ml maximize, difficult

nlcom exp([lndiff]_cons)+[cut1]_cons

ml model lf myoprobit2_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_pminorty_dw_z, nocons) (cut1:) (lndiff:) /*
*/(xb2:f_c9010t_pminorty_dw_z = x_f_release1 $x_covariates interaction*) (lnsigma:) (rho:) [pw=$wt] , /*
*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
*/diparm(lnsigma, exp label(sigma)) 
ml maximize

nlcom exp([lndiff]_cons)+[cut1]_cons

ml model lf myoprobit3_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z, nocons) (cut1:) (lndiff:) /*
						*/(xb2:f_c9010t_perpov_dw_z = x_f_release1 $x_covariates interaction*) (lnvar_v1:) (rho_uv1:) /*
						*/(xb3:f_c9010t_pminorty_dw_z = x_f_release1 $x_covariates interaction*) (lnvar_v2:) (rho_uv2:) (cov_v1v2:) [pw=$wt], /*
						*/diparm(lndiff cut1, f(exp(@1)+@2) d(exp(@1) -exp(@1)+1) label(cut2)) /*
						*/diparm(lnvar_v1, exp label(var_v1)) /*
						*/diparm(lnvar_v2, exp label(var_v2)) 
ml init xb1:x_f_site_bal=0.3039 x_f_site_bos=0.3386 x_f_site_chi=0.1527 x_f_site_la=0.0639 f_c9010t_perpov_dw_z=-0.2998 f_c9010t_pminorty_dw_z=0.3151 cut1:_cons=-0.8165
ml maximize, difficult

nlcom exp([lndiff]_cons)+[cut1]_cons
nlcom (exp([lnvar_v1]_cons))^0.5
nlcom (exp([lnvar_v2]_cons))^0.5
