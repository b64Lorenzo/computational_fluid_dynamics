function [coeffs,offsets,source] = b3_scheme_1d(mesh,U,k,i)

[coeffs,offsets,source] = lambda_scheme_1d( mesh,U,k,i,0.5);

end