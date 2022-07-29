
** TABLE D1

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
using "${pathres}tableD1_${date}.tex", cells("b(star label() fmt(4))" "se(par label() fmt(4))") ///
replace mlabel(none) collabels(none) eqlabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) stardetach label ///
varwidth(12) modelwidth(10) style(tex) drop($site_covs x_f_release1 ${x_covariates}  f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z cut1:_cons lndiff:_cons rho:_cons xb2:x_f_release1 xb2:$x_covariates xb2:interaction_3_5 lnvar_v1:_cons rho_uv1:_cons lnvar_v2:_cons rho_uv2:_cons xb3:x_f_release1 xb3:$x_covariates xb3:interaction_3_5 lnsigma:_cons cov_v1v2:_cons) stats(N, fmt(%9.0f)) ///
order(xb2:interaction_1_1 xb2:interaction_1_2 xb2:interaction_1_3 xb2:interaction_1_4 xb2:interaction_1_5 xb2:interaction_2_1 xb2:interaction_2_2 xb2:interaction_2_3 xb2:interaction_2_4 xb2:interaction_2_5 xb2:interaction_3_1 xb2:interaction_3_2 xb2:interaction_3_3 xb2:interaction_3_4 xb2:_cons ///
xb3:interaction_1_1 xb3:interaction_1_2 xb3:interaction_1_3 xb3:interaction_1_4 xb3:interaction_1_5 xb3:interaction_2_1 xb3:interaction_2_2 xb3:interaction_2_3 xb3:interaction_2_4 xb3:interaction_2_5 xb3:interaction_3_1 xb3:interaction_3_2 xb3:interaction_3_3 xb3:interaction_3_4 xb3:_cons) ///
varlabels(xb2:interaction_1_1 "$ \delta^{exp, Balt}_{1} $" xb2:interaction_1_2 "$ \delta^{exp, Bos}_{1} $" xb2:interaction_1_3 "$ \delta^{exp, Chi}_{1} $" xb2:interaction_1_4 "$ \delta^{exp, LA}_{1} $" xb2:interaction_1_5 "$ \delta^{exp, NY}_{1} $" ///
xb2:interaction_2_1 "$ \delta^{sec8, Balt}_{1} $" xb2:interaction_2_2 "$ \delta^{sec8, Bos}_{1} $" xb2:interaction_2_3 "$ \delta^{sec8, Chi}_{1} $" xb2:interaction_2_4 "$ \delta^{sec8, LA}_{1} $" xb2:interaction_2_5 "$ \delta^{sec8, NY}_{1} $" ///
xb2:interaction_3_1 "$ \delta^{cont, Balt}_{1} $" xb2:interaction_3_2 "$ \delta^{cont, Bos}_{1} $" xb2:interaction_3_3 "$ \delta^{cont, Chi}_{1} $" xb2:interaction_3_4 "$ \delta^{cont, LA}_{1} $" xb2:_cons "$ \delta^{cont, NY}_{1} $" ///
xb3:interaction_1_1 "$ \delta^{exp, Balt}_{2} $" xb3:interaction_1_2 "$ \delta^{exp, Bos}_{2} $" xb3:interaction_1_3 "$ \delta^{exp, Chi}_{2} $" xb3:interaction_1_4 "$ \delta^{exp, LA}_{2} $" xb3:interaction_1_5 "$ \delta^{exp, NY}_{2} $" ///
xb3:interaction_2_1 "$ \delta^{sec8, Balt}_{2} $" xb3:interaction_2_2 "$ \delta^{sec8, Bos}_{2} $" xb3:interaction_2_3 "$ \delta^{sec8, Chi}_{2} $" xb3:interaction_2_4 "$ \delta^{sec8, LA}_{2} $" xb3:interaction_2_5 "$ \delta^{sec8, NY}_{2} $" ///
xb3:interaction_3_1 "$ \delta^{cont, Balt}_{2} $" xb3:interaction_3_2 "$ \delta^{cont, Bos}_{2} $" xb3:interaction_3_3 "$ \delta^{cont, Chi}_{2} $" xb3:interaction_3_4 "$ \delta^{cont, LA}_{2} $" xb3:_cons "$ \delta^{cont, NY}_{2} $") ///
prehead("\begin{ThreePartTable}\begin{center} \caption{@title}\renewcommand{\arraystretch}{1.2} \small " ///
"\begin{TableNotes} \scriptsize " ///
"\item \textit{Notes:} Each column reports first stage estimates of a triangular model for specifications reported in corresponding columns of Table \ref{tab:tiv}. The dependent variables in the first stage are neighborhood poverty and neighborhood minority. Columns (1)-(3) exclude while columns (4)-(6) include a complete set of baseline characteristics (as given in Table \ref{tab:table1}), as well as whether a sample adult was included in the first release of the long-term evaluation survey fielding period, as covariates X; all regressions are weighted; * p-value $<0.10$, ** p-value $<0.05$, *** p-value $<0.01$." ///
"\item \textit{Source:} Data from ICPSR Study 34860: Moving to Opportunity: Final Impacts Evaluation Science Article Data, 2008-2010." ///
"\end{TableNotes}" /// 
"\begin{longtable}{@{}l r@{}l r@{}l r@{}l r@{}l r@{}l r@{}l} " ///
"\caption{Triangular IV estimation of neighborhood effects on SWB, first stage} \\ \toprule \label{tab:tiv_firststage} " ///
" & \multicolumn{2}{c}{(1)} & \multicolumn{2}{c}{(2)} & \multicolumn{2}{c}{(3)} & \multicolumn{2}{c}{(4)} & \multicolumn{2}{c}{(5)} & \multicolumn{2}{c}{(6)}\\") ///
posthead("\midrule \endhead \hline \multicolumn{13}{r}{\textit{continued on next page}} \endfoot \endlastfoot  ") prefoot("\midrule ") ///
postfoot("\bottomrule \end{longtable}" ///
" \insertTableNotes" ///
"\end{center} \end{ThreePartTable}")
