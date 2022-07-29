
*codebook happy_scale012 f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z ${site_covs} ${x_covariates}

regress happy_scale012 f_c9010t_perpov_dw_z f_c9010t_pminorty_dw_z ${site_covs} ${x_covariates}
gen sample_allx = e(sample)

sort ra_group

by ra_group: eststo: estpost su ${site_covs} sgx_rasite_3g_all_nyc if sample_allx==1, listwise

esttab using ${pathres}tableB1_${date}.tex, replace main(mean) unstack noobs cells("mean(fmt(2))") ///
label booktabs nonum f collabels(none) gaps plain ///
varlabels(x_f_site_balt "Baltimore" x_f_site_bos "Boston" x_f_site_chi "Chicago" x_f_site_la "Los Angeles" sgx_rasite_3g_all_nyc "New York") ///
prehead("\begin{ThreePartTable}\begin{center}\caption{@title}\renewcommand{\arraystretch}{1.2}\small" ///
"\begin{TableNotes} \scriptsize " ///
"\item \textit{Notes:} Each cell gives the average value of a variable in the sub-sample. Only observations with non-missing values for Subjective Well Being (SWB), neighbourhood characteristics and x covariates are used. There are 7/3,273 observations with missing SWB, 3/3,273 observations with missing neighborhood characteristics and 89/3,273 observations with missing household income. Some observations of x covariates include imputed values." ///
"\item \textit{Source:} Data from ICPSR Study 34860: Moving to Opportunity: Final Impacts Evaluation Science Article Data, 2008-2010." ///
"\end{TableNotes}" /// 
"\begin{longtable}{l c c c c} " ///
"\caption{Baseline characteristics of MTO adults or covariates X across randomization groups} \\ \toprule \label{tab:table1} ") ///
posthead("\midrule \endhead \hline \multicolumn{5}{r}{\textit{continued on next page}} \endfoot \endlastfoot  \textit{\textbf{Site:}} & & & & \\")  

eststo clear

by ra_group: eststo: estpost su x_rad_ad_ethrace_black_nh x_rad_ad_ethrace_hisp tmp_female x_rad_ad_le_35 x_rad_ad_36_40 x_rad_ad_41_45 x_rad_ad_46_50 ///
x_f_ad_nevmarr x_f_ad_parentu18 x_f_ad_working x_f_ad_edinsch x_f_ad_edgradhs x_f_ad_edged x_f_hh_afdc if sample_allx==1, listwise

esttab using ${pathres}tableB1_${date}.tex, append main(mean) unstack noobs nomtitles cells("mean(fmt(2))") ///
label booktabs nonum f collabels(none) gaps plain ///
varlabels(x_rad_ad_ethrace_black_nh "African American (non-hispanic)" x_rad_ad_ethrace_hisp "Hispanic ethnicity (any race)" tmp_female "Female" ///
x_rad_ad_le_35 "$<=35$ years old" x_rad_ad_36_40 "36-40 years old" x_rad_ad_41_45 "41-45 years old" x_rad_ad_46_50 "46-50 years old" ///
x_f_ad_nevmarr "Never married" x_f_ad_parentu18 "Parent while younger than 18 years old" x_f_ad_working "Working" x_f_ad_edinsch "Enrolled in school" ///
x_f_ad_edgradhs "High school diploma" x_f_ad_edged "General Education Development (GED) certificate" ///
x_f_hh_afdc "Receiving Aid to Families with Dependent Children (AFDC)") /// 
prehead(" \\ \textit{\textbf{Demographic characteristics:}} & & & & \\")

eststo clear

by ra_group: eststo: estpost su rad_svy_bl_totincm_2009d if sample_allx==1, listwise

esttab using ${pathres}tableB1_${date}.tex, append main(mean) unstack noobs nomtitles cells("mean(fmt(%8.0gc))") ///
label booktabs nonum f collabels(none) gaps plain ///
varlabels(rad_svy_bl_totincm_2009d "Household income (dollars)") ///
prehead(" \\ \textit{\textbf{Household characteristics:}} & & & & \\")

eststo clear

by ra_group: eststo: estpost su x_f_hh_car x_f_hh_disabl x_f_hh_noteens x_f_hh_size2 x_f_hh_size3 x_f_hh_size4 if sample_allx==1, listwise

esttab using ${pathres}tableB1_${date}.tex, append main(mean) unstack noobs nomtitles cells("mean(fmt(2))") ///
label booktabs nonum f collabels(none) gaps plain ///
varlabels(x_f_hh_car "Household owns a car" x_f_hh_disabl "Household member had a disability" x_f_hh_noteens "No teens in household" ///
x_f_hh_size2 "Household size is $<=2$" x_f_hh_size3 "Household size is 3" x_f_hh_size4 "Household size is 4") 

eststo clear

by ra_group: eststo: estpost su x_f_hh_victim x_f_hood_unsafenit x_f_hood_verydissat x_f_hood_5y x_f_hous_mov3tm  ///
x_f_hood_nofamily x_f_hood_nofriend x_f_hood_chat x_f_hood_nbrkid x_f_hous_fndapt x_f_hous_sec8bef if sample_allx==1, listwise

esttab using ${pathres}tableB1_${date}.tex, append main(mean) unstack noobs nomtitles cells("mean(fmt(2))") ///
label booktabs nonum f collabels(none) gaps plain ///
varlabels(x_f_hh_victim "Household member was a crime victim in past 6 months" x_f_hood_unsafenit "Neighborhood streets very unsafe at night" ///
x_f_hood_verydissat "Very dissatisfied with neighborhood" x_f_hood_5y "Household living in neighborhood $>5$ years" ///
x_f_hous_mov3tm "Household moved more $>3x$ in last 5 yrs" x_f_hood_nofamily "Household has no family living in neighborhood" ///
x_f_hood_nofriend "Household has no friends living in neighborhood" x_f_hood_chat "Household head chatted with neighbor $>=1x$ per week" ///
x_f_hood_nbrkid "Household head very likely to to tell on neighborhood kid" x_f_hous_fndapt "Household head very sure of finding apartment" ///
x_f_hous_sec8bef "Housheold head applied for Section 8 before") ///
prehead(" \\ \textit{\textbf{Neighborhood characteristics:}} & & & & \\")

eststo clear

by ra_group: eststo: estpost su x_f_hous_movdrgs x_f_hous_movschl cov_hous_movapt cov_hous_movjob if sample_allx==1, listwise

esttab using ${pathres}tableB1_${date}.tex, append main(mean) unstack nomtitles cells("mean(fmt(2))") ///
label booktabs nonum f collabels(none) gaps plain stats(N, fmt(%9.0f)) ///
varlabels(x_f_hous_movdrgs "Want to move to get away from gangs and drugs" x_f_hous_movschl "Want to move for better schools for children" ///
cov_hous_movapt "Want to move to get a bigger/better apartment" cov_hous_movjob "Want to move to get a job") ///
prehead(" \\ \textit{\textbf{Primary or secondary reason for wanting to move:}} & & & & \\") ///
prefoot("\midrule ") postfoot( ///
"\bottomrule \end{longtable}" ///
" \insertTableNotes" ///
" \end{center} \end{ThreePartTable}")
