function [a_W,a_P,a_E,source] = central_scheme_1d(mesh,U,k,i)

%==========================================================================
%
% CENTRAL_SCHEME
%
% Central differencing discretization of
%
%     U d(phi)/dx = k d²(phi)/dx²
%
%==========================================================================
%
% INPUT
%
% mesh.x  : grid coordinates
% U       : convection velocity
% k       : diffusion coefficient
% i       : node index
%
% OUTPUT
%
% a_W     : west coefficient
% a_P     : central coefficient
% a_E     : east coefficient
% source  : source term contribution
%
%==========================================================================

x = mesh.x;

dx_W = x(i)   - x(i-1);
dx_E = x(i+1) - x(i);

D_W = k / dx_W;
D_E = k / dx_E;

F = U;

a_W = D_W + F/2;

a_E = D_E - F/2;

a_P = -(a_W + a_E);

source = 0;

end