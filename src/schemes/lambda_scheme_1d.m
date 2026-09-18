function [coeffs,offsets,source] = ...
    lambda_scheme_1d(mesh,U,k,i,lambda)

x = mesh.x(:);

h = x(2)-x(1);

if max(abs(diff(x)-h)) > 1e-12
    error(['Lambda schemes require uniform grids']);
end

c_im2 = U*(lambda/h);

c_im1 = -U*(6*lambda+1)/(2*h) -k/h^2;

c_i   = U*(1+6*lambda)/(2*h)+2*k/h^2;

c_ip1 = -U*(lambda/h) -k/h^2;

coeffs  = [c_im2 c_im1 c_i c_ip1];
offsets = [-2 -1 0 1];
source  = 0;

end