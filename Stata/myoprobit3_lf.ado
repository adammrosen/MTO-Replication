program myoprobit3_lf
			args lnf xb1 a1 lndiff xb2 lnvar_v1 rho_uv1 xb3 lnvar_v2 rho_uv2 cov_v1v2
			
			tempvar lambdav1 lambdav2 cov_uv1 cov_uv2 detsigmav mu1 mu2 sigma_p sigma
			
			qui{
			gen double `lambdav1'	= $ML_y2 - `xb2'
			gen double `lambdav2'	= $ML_y3 - `xb3'
			
			gen double `cov_uv1'	= `rho_uv1'*((exp(`lnvar_v1'))^0.5)
			gen double `cov_uv2'	= `rho_uv2'*((exp(`lnvar_v2'))^0.5)

			gen double `detsigmav'	= (exp(`lnvar_v1')*exp(`lnvar_v2'))-((`cov_v1v2')^2)

			gen double `mu1' 		= (((`cov_uv1')*exp(`lnvar_v2'))-((`cov_uv2')*(`cov_v1v2')))/`detsigmav'
			gen double `mu2' 		= (((`cov_uv1')*(`cov_v1v2'))-((`cov_uv2')*exp(`lnvar_v1')))/`detsigmav'

			gen double `sigma_p'	= ((((`cov_uv1'))^2)*exp(`lnvar_v2')) - (((2*(`cov_uv1'))*(`cov_uv2'))*(`cov_v1v2')) + ((((`cov_uv2'))^2)*exp(`lnvar_v1'))
			gen double `sigma' 		= (1-(`sigma_p'/`detsigmav'))^0.5
			
			replace `lnf' = ln(normal((`a1'-`xb1'-((`mu1'*`lambdav1')-(`mu2'*`lambdav2')))/`sigma')) + ln((1/((2*3.1415926)*((`detsigmav')^0.5)))*exp(-0.5* (((((`lambdav1')^2)*exp(`lnvar_v2'))-2*((`lambdav1'*`lambdav2')*(`cov_v1v2'))+(((`lambdav2')^2)*exp(`lnvar_v1')))/(`detsigmav')))) if $ML_y1 == 0
			replace `lnf' = ln(normal((exp(`lndiff')+`a1'-`xb1'-((`mu1'*`lambdav1')-(`mu2'*`lambdav2')))/`sigma')-normal((`a1'-`xb1'-((`mu1'*`lambdav1')-(`mu2'*`lambdav2')))/`sigma'))+ ln((1/((2*3.1415926)*((`detsigmav')^0.5)))*exp(-0.5* (((((`lambdav1')^2)*exp(`lnvar_v2'))-2*((`lambdav1'*`lambdav2')*(`cov_v1v2'))+(((`lambdav2')^2)*exp(`lnvar_v1')))/(`detsigmav')))) if $ML_y1 == 1
			replace `lnf' = ln(1-normal((exp(`lndiff')+`a1'-`xb1'-((`mu1'*`lambdav1')-(`mu2'*`lambdav2')))/`sigma'))+ ln((1/((2*3.1415926)*((`detsigmav')^0.5)))*exp(-0.5* (((((`lambdav1')^2)*exp(`lnvar_v2'))-2*((`lambdav1'*`lambdav2')*(`cov_v1v2'))+(((`lambdav2')^2)*exp(`lnvar_v1')))/(`detsigmav')))) if $ML_y1 == 2
			
			}
end
