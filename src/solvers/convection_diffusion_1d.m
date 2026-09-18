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
%     - Classical 3-point schemes (TDMA)
%     - Lambda/QUICK/B3 5-diagonal schemes (PDMA)
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

n_nodes    = length(x);
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

a_WW = zeros(n_internal,1);  % i-2
a_W  = zeros(n_internal,1);  % i-1
a_P  = zeros(n_internal,1);  % i
a_E  = zeros(n_internal,1);  % i+1
a_EE = zeros(n_internal,1);  % i+2

rhs = zeros(n_internal,1);

has_pentadiagonal = false;

%--------------------------------------------------------------------------
% Assembly loop
%--------------------------------------------------------------------------

for I = 1:n_internal

    i = I + 1;

    try

        %--------------------------------------------------------------
        % New lambda/QUICK/B3 format
        %
        % [coeffs,offsets,source]
        %--------------------------------------------------------------

        [coeffs,offsets,source] = ...
            scheme(mesh,U,k,i);

        if ~isempty(offsets)

            for m = 1:length(offsets)

                switch offsets(m)

                    case -2
                        a_WW(I) = coeffs(m);

                    case -1
                        a_W(I) = coeffs(m);

                    case 0
                        a_P(I) = coeffs(m);

                    case 1
                        a_E(I) = coeffs(m);

                    case 2
                        a_EE(I) = coeffs(m);

                    otherwise
                        error('Unsupported stencil offset.');
                end

            end

            if any(abs(offsets) == 2)
                has_pentadiagonal = true;
            end

        end

    catch

        %--------------------------------------------------------------
        % Old format
        %
        % [aW,aP,aE,source]
        %--------------------------------------------------------------

        [coeff_W,coeff_P,coeff_E,source] = ...
            scheme(mesh,U,k,i);

        a_W(I) = coeff_W;
        a_P(I) = coeff_P;
        a_E(I) = coeff_E;

    end

    rhs(I) = source;

end

%--------------------------------------------------------------------------
% West boundary condition
%--------------------------------------------------------------------------

if strcmpi(bc.W.type,'dirichlet')

    %
    % First interior equation
    %

    rhs(1) = rhs(1) ...
           - a_W(1)*bc.W.value;

    %
    % Second interior equation may contain phi(1)
    % through the WW coefficient
    %

    if has_pentadiagonal && n_internal >= 2

        rhs(2) = rhs(2) ...
               - a_WW(2)*bc.W.value;

    end

elseif strcmpi(bc.W.type,'neumann')

    dx_W = x(2) - x(1);

    %
    % Original tridiagonal treatment
    %

    a_P(1) = a_P(1) + a_W(1);

    rhs(1) = rhs(1) ...
           + a_W(1)*bc.W.value*dx_W;

    if has_pentadiagonal

        warning(['Pentadiagonal scheme with ', ...
                 'Neumann BC not fully generalized.']);

    end

else

    error('Unknown west boundary condition type.');

end

%--------------------------------------------------------------------------
% East boundary condition
%--------------------------------------------------------------------------

if strcmpi(bc.E.type,'dirichlet')

    %
    % Last interior equation
    %

    rhs(end) = rhs(end) ...
             - a_E(end)*bc.E.value;

    %
    % Penultimate equation may contain phi(end)
    %

    if has_pentadiagonal && n_internal >= 2

        rhs(end-1) = rhs(end-1) ...
                   - a_EE(end-1)*bc.E.value;

    end

elseif strcmpi(bc.E.type,'neumann')

    dx_E = x(end) - x(end-1);

    a_P(end) = a_P(end) + a_E(end);

    rhs(end) = rhs(end) ...
             - a_E(end)*bc.E.value*dx_E;

    if has_pentadiagonal

        warning(['Pentadiagonal scheme with ', ...
                 'Neumann BC not fully generalized.']);

    end

else

    error('Unknown east boundary condition type.');

end

%--------------------------------------------------------------------------
% Solve system
%--------------------------------------------------------------------------

if has_pentadiagonal

    logger(fid,'INFO', ...
        'Pentadiagonal stencil detected -> using PDMA');

    %
    % PDMA storage
    %
    % b2 = diag(A,-2)
    % b1 = diag(A,-1)
    % a  = diag(A,0)
    % c1 = diag(A,+1)
    % c2 = diag(A,+2)
    %

    b2 = zeros(n_internal,1);
    b1 = zeros(n_internal,1);
    a  = a_P;
    c1 = zeros(n_internal,1);
    c2 = zeros(n_internal,1);

    b2(3:end) = a_WW(3:end);
    b1(2:end) = a_W(2:end);

    c1(1:end-1) = a_E(1:end-1);
    c2(1:end-2) = a_EE(1:end-2);

    %
    % Store matrix diagonals
    %

    A.a_WW = a_WW;
    A.a_W  = a_W;
    A.a_P  = a_P;
    A.a_E  = a_E;
    A.a_EE = a_EE;

    A.b2 = b2;
    A.b1 = b1;
    A.a  = a;
    A.c1 = c1;
    A.c2 = c2;

    phi_internal = pdma( ...
        b2,...
        b1,...
        a,...
        c1,...
        c2,...
        rhs,...
        n_internal);

else

    logger(fid,'INFO', ...
        'Tridiagonal stencil detected -> using TDMA');

    lower_diag = zeros(n_internal,1);
    main_diag  = a_P;
    upper_diag = zeros(n_internal,1);

    lower_diag(2:end)  = a_W(2:end);
    upper_diag(1:end-1)= a_E(1:end-1);

    A.a_W = a_W;
    A.a_P = a_P;
    A.a_E = a_E;

    A.lower_diag = lower_diag;
    A.main_diag  = main_diag;
    A.upper_diag = upper_diag;

    phi_internal = tdma( ...
        lower_diag,...
        main_diag,...
        upper_diag,...
        rhs,...
        n_internal);

end

%--------------------------------------------------------------------------
% Assemble full solution
%--------------------------------------------------------------------------

phi = zeros(n_nodes,1);

phi(2:end-1) = phi_internal;

if strcmpi(bc.W.type,'dirichlet')

    phi(1) = bc.W.value;

else

    dx_W = x(2)-x(1);

    phi(1) = phi(2) ...
           - bc.W.value*dx_W;

end

if strcmpi(bc.E.type,'dirichlet')

    phi(end) = bc.E.value;

else

    dx_E = x(end)-x(end-1);

    phi(end) = phi(end-1) ...
             + bc.E.value*dx_E;

end

logger(fid,'INFO', ...
    'Solution assembled (%d total nodes)', ...
    n_nodes);

logger(fid,'INFO', ...
    'Convection-diffusion solution completed');

end