
** TABLE 1

ml model lf myoprobit_lf (xb1:happy_scale012 = $site_covs f_c9010t_perpov_dw_z, nocons) (cut1:) (cut2:) [pw=$wt]
ml maximize
est store pov_op_nox

ml model lf myoprobit_lf (xb1:happy_scale012 = $site_covs f_c9010t_pminorty_dw_z, nocons) (cut1:) (cut2:) [pw=$wt]
ml maximize
est store min_op_nox

ml model lf myoprobit_lf (xb1:happy_scale012 = $site_covs f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z, nocons) (cut1:) (cut2:) [pw=$wt]
ml maximize
est store povmin_op_nox

ml model lf myoprobit_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z, nocons) (cut1:) (cut2:) [pw=$wt]
ml maximize
est store pov_op_x

ml model lf myoprobit_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_pminorty_dw_z, nocons) (cut1:) (cut2:) [pw=$wt]
ml maximize
est store min_op_x

ml model lf myoprobit_lf (xb1:happy_scale012 = $site_covs x_f_release1 $x_covariates f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z, nocons) (cut1:) (cut2:) [pw=$wt]
ml maximize
est store povmin_op_x

estout pov_op_nox min_op_nox povmin_op_nox pov_op_x min_op_x povmin_op_x ///
using "${pathres}table1_${date}.tex", cells("b(star label() fmt(4))" "se(par label() fmt(4))") ///
replace mlabel(none) collabels(none) eqlabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) stardetach label ///
varwidth(12) modelwidth(10) style(tex) drop(${x_covariates} x_f_release1) stats(N, fmt(%9.0f)) ///
order(f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z x_f_site_balt x_f_site_bos x_f_site_chi x_f_site_la) ///
varlabels(f_c9010t_perpov_dw_z "$ \beta_{Poverty} $" f_c9010t_pminorty_dw_z "$ \beta_{Minority} $" x_f_site_balt "$ \gamma_{Baltomore} $" x_f_site_bos "$ \gamma_{Boston} $" ///
x_f_site_chi "$ \gamma_{Chicago} $" x_f_site_la "$ \gamma_{LA} $" cut1:_cons "$ c_{1} $" cut2:_cons "$ c_{2} $") ///
title(Ordered probit estimation of neighborhood effects on SWB \label{tab:op}) ///
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
