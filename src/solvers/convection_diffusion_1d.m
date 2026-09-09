function [phi,A,rhs] = convection_diffusion_1d( ...
    mesh,U,k,bc,scheme,fid)

%==========================================================================
%
% CONVECTION_DIFFUSION_1D
%
% Solves:
%
%     U d(phi)/dx = k d²(phi)/dx²
%
% on a uniform or non-uniform grid.
%
% Supports:
%
%     - Dirichlet BCs
%     - Neumann BCs
%
%==========================================================================

if nargin < 6
    fid = [];
end

if nargin < 4 || isempty(bc)

    bc.W.type  = 'dirichlet';
    bc.W.value = 0;

    bc.E.type  = 'dirichlet';
    bc.E.value = 0;

end

x = mesh.x(:);

n_nodes = length(x);

n_internal = n_nodes - 2;

if n_internal < 1
    error('Mesh must contain at least 3 nodes.');
end

logger(fid,'INFO', ...
    'Discretization scheme: %s', ...
    func2str(scheme));

logger(fid,'INFO', ...
    'Assembling convection-diffusion system (%d internal nodes)', ...
    n_internal);

%--------------------------------------------------------------------------
% Allocate storage
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
% West boundary condition
%--------------------------------------------------------------------------

if strcmpi(bc.W.type,'dirichlet')

    rhs(1) = rhs(1) - a_W(1)*bc.W.value;

elseif strcmpi(bc.W.type,'neumann')

    dx_W = x(2) - x(1);

    a_P(1) = a_P(1) + a_W(1);

    rhs(1) = rhs(1) + a_W(1)*bc.W.value*dx_W;

else

    error('Unknown west boundary condition type.');

end

%--------------------------------------------------------------------------
% East boundary condition
%--------------------------------------------------------------------------

if strcmpi(bc.E.type,'dirichlet')

    rhs(end) = rhs(end) - a_E(end)*bc.E.value;

elseif strcmpi(bc.E.type,'neumann')

    dx_E = x(end) - x(end-1);

    a_P(end) = a_P(end) + a_E(end);

    rhs(end) = rhs(end) - a_E(end)*bc.E.value*dx_E;

else

    error('Unknown east boundary condition type.');

end

%--------------------------------------------------------------------------
% Construct and storing TDMA diagonals and Solving
%--------------------------------------------------------------------------

lower_diag = zeros(n_internal,1);
main_diag  = a_P;
upper_diag = zeros(n_internal,1);

lower_diag(2:end) = a_W(2:end);

upper_diag(1:end-1) = a_E(1:end-1);

A.a_W = a_W;
A.a_P = a_P;
A.a_E = a_E;

A.lower_diag = lower_diag;
A.main_diag  = main_diag;
A.upper_diag = upper_diag;

logger(fid,'INFO','Solving tridiagonal system using TDMA');

phi_internal = tdma( ...
    lower_diag,...
    main_diag,...
    upper_diag,...
    rhs,...
    n_internal);

phi = zeros(n_nodes,1);

phi(2:end-1) = phi_internal;

if strcmpi(bc.W.type,'dirichlet')

    phi(1) = bc.W.value;

else

    dx_W = x(2) - x(1);

    phi(1) = phi(2) - bc.W.value*dx_W;

end

if strcmpi(bc.E.type,'dirichlet')

    phi(end) = bc.E.value;

else

    dx_E = x(end) - x(end-1);

    phi(end) = phi(end-1) + bc.E.value*dx_E;

end

logger(fid,'INFO', ...
    'Solution assembled (%d total nodes)', ...
    n_nodes);

logger(fid,'INFO', ...
    'Convection-diffusion solution completed');

end