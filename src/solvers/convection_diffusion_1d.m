function [phi,A,rhs] = convection_diffusion_1d( ...
    mesh,U,k,phi_W,phi_E,scheme,fid)

%==========================================================================
%
% CONVECTION_DIFFUSION_1D
%
% Solves the steady one-dimensional convection-diffusion equation
%
%     U d(phi)/dx = k d²(phi)/dx²
%
% on a uniform or non-uniform grid.
%
% Boundary conditions:
%
%     phi(0) = phi_W
%     phi(L) = phi_E
%
% INPUT
%
% mesh.x      : grid coordinates
% U           : convection velocity
% k           : diffusion coefficient
% phi_W       : west boundary value
% phi_E       : east boundary value
% scheme      : discretization scheme function handle
% fid         : logger file identifier
%
% OUTPUT
%
% phi         : solution including boundary nodes
% A           : structure containing matrix diagonals
% rhs         : right-hand-side vector
%
% Example:
%
% [phi,A,rhs] = convection_diffusion_1d( ...
%     mesh,U,k,phi_W,phi_E,@central_scheme,fid);
%
%==========================================================================

if nargin < 7
    fid = [];
end

x = mesh.x(:);

n_nodes = length(x);

n_internal = n_nodes - 2;

if n_internal < 1
    error('Mesh must contain at least 3 nodes.');
end

logger(fid,'INFO', ...
    'Assembling convection-diffusion system (%d internal nodes)', ...
    n_internal);

%--------------------------------------------------------------------------
% Tridiagonal matrix coefficients
%--------------------------------------------------------------------------

a_W = zeros(n_internal,1);
a_P = zeros(n_internal,1);
a_E = zeros(n_internal,1);

rhs = zeros(n_internal,1);

%--------------------------------------------------------------------------
% Assembly loop
%--------------------------------------------------------------------------

for I = 1:n_internal

    i = I + 1;

    [coeff_W,coeff_P,coeff_E,source] = ...
        scheme(mesh,U,k,i);

    a_W(I) = coeff_W;
    a_P(I) = coeff_P;
    a_E(I) = coeff_E;

    rhs(I) = source;

end

%--------------------------------------------------------------------------
% Apply Dirichlet boundary conditions
%--------------------------------------------------------------------------

rhs(1) = rhs(1) - a_W(1)*phi_W;

rhs(end) = rhs(end) - a_E(end)*phi_E;

%--------------------------------------------------------------------------
% Construct TDMA diagonals
%--------------------------------------------------------------------------

lower = zeros(n_internal,1);
diag  = a_P;
upper = zeros(n_internal,1);

lower(2:end) = a_W(2:end);

upper(1:end-1) = a_E(1:end-1);

A.lower = lower;
A.diag  = diag;
A.upper = upper;

logger(fid,'INFO','Solving tridiagonal system using TDMA');

%--------------------------------------------------------------------------
% Solve system
%--------------------------------------------------------------------------

phi_internal = tdma(lower,diag,upper,rhs,n_internal);

%--------------------------------------------------------------------------
% Assemble full solution vector
%--------------------------------------------------------------------------

phi = zeros(n_nodes,1);

phi(1) = phi_W;

phi(end) = phi_E;

phi(2:end-1) = phi_internal;

logger(fid,'INFO', ...
    'Convection-diffusion solution completed');

end