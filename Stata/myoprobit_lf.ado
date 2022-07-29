program myoprobit_lf
			args lnfj xb a1 a2
			
			quietly replace `lnfj' = ln(  normal(`a1'-`xb')) if $ML_y1 == 0
			quietly replace `lnfj' = ln(  normal(`a2'-`xb') - normal(`a1'-`xb')) if $ML_y1 == 1
			quietly replace `lnfj' = ln(1-normal(`a2'-`xb')) if $ML_y1 == 2
end
