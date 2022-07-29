
** FIGURE 1

kdensity f_c9010t_perpov_dw_z, kernel(gaussian) nograph generate(x1 fx1)
kdensity f_c9010t_perpov_dw_z if ra_group==1, kernel(gaussian) nograph generate(fx1_1) at(x1)
kdensity f_c9010t_perpov_dw_z if ra_group==2, kernel(gaussian) nograph generate(fx1_2) at(x1)
kdensity f_c9010t_perpov_dw_z if ra_group==3, kernel(gaussian) nograph generate(fx1_3) at(x1)
label var fx1_1 "Experimental"
label var fx1_2 "Section 8"
label var fx1_3 "Control"
line fx1_1 fx1_2 fx1_3 x1, xtitle("Neighborhood Poverty, z-score") ytitle("Density") lpattern(dash dot solid) lcolor(black black black) ylab( , grid) scale(0.8)
graph export ${pathres}w1_density_${date}.eps, replace

count if ra_group==1&missing(f_c9010t_perpov_dw_z)==0
count if ra_group==2&missing(f_c9010t_perpov_dw_z)==0
count if ra_group==3&missing(f_c9010t_perpov_dw_z)==0
count if missing(f_c9010t_perpov_dw_z)

** FIGURE 2

kdensity f_c9010t_pminorty_dw_z, kernel(gaussian) nograph generate(x2 fx2)
kdensity f_c9010t_pminorty_dw_z if ra_group==1, kernel(gaussian) nograph generate(fx2_1) at(x2)
kdensity f_c9010t_pminorty_dw_z if ra_group==2, kernel(gaussian) nograph generate(fx2_2) at(x2)
kdensity f_c9010t_pminorty_dw_z if ra_group==3, kernel(gaussian) nograph generate(fx2_3) at(x2)
label var fx2_1 "Experimental"
label var fx2_2 "Section 8"
label var fx2_3 "Control"
line fx2_1 fx2_2 fx2_3 x2, xtitle("Neighborhood Minority, z-score") ytitle("Density") lpattern(dash dot solid) lcolor(black black black) ylab( , grid) scale(0.8)
graph export ${pathres}w2_density_${date}.eps, replace

count if ra_group==1&missing(f_c9010t_pminorty_dw_z)==0
count if ra_group==2&missing(f_c9010t_pminorty_dw_z)==0
count if ra_group==3&missing(f_c9010t_pminorty_dw_z)==0
count if missing(f_c9010t_pminorty_dw_z)
