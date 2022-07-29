program myoprobit2_lf
			args lnf xb1 a1 lndiff xb2 lnsigma rho
			
			tempvar lambda N D
			
			qui{
			gen double `lambda'	= $ML_y2 - `xb2'
			
			gen double `N' 		= (`rho')/exp(`lnsigma')
			gen double `D' 		= (1-(`rho')^2)^0.5
			
			replace `lnf' = ln(normal((`a1'-`xb1'-(`N'*`lambda'))/`D')) + ln(normalden($ML_y2 ,`xb2',exp(`lnsigma'))) if $ML_y1 == 0
			replace `lnf' = ln(normal((exp(`lndiff')+`a1'-`xb1'-(`N'*`lambda'))/`D') - normal((`a1'-`xb1'-(`N'*`lambda'))/`D')) + ln(normalden($ML_y2 ,`xb2',exp(`lnsigma'))) if $ML_y1 == 1
			replace `lnf' = ln(1-normal((exp(`lndiff')+`a1'-`xb1'-(`N'*`lambda'))/`D')) + ln(normalden($ML_y2 ,`xb2',exp(`lnsigma'))) if $ML_y1 == 2
			}
end
