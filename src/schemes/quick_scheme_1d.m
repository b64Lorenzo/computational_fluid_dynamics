function [coeffs,offsets,source] = quick_scheme_1d(mesh,U,k,i)

[coeffs,offsets,source] = lambda_scheme_1d(mesh,U,k,i,1/8);

end