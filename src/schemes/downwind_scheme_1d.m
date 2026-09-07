function [a_im1,a_i,a_ip1,source] = downwind_scheme_1d(mesh,U,k,i)

x = mesh.x;

dx_im1_i = x(i)   - x(i-1);
dx_i_ip1 = x(i+1) - x(i);

D_im1 = k/dx_im1_i;
D_ip1 = k/dx_i_ip1;

a_im1 = D_im1 + min(U,0);

a_ip1 = D_ip1 + max(U,0);

a_i = -(a_im1 + a_ip1);

source = 0;

end