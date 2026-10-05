function [A,rhs,has_pentadiagonal] = ...
    assemble_convection_diffusion_1d( ...
    mesh,U,k,bc,scheme)

x = mesh.x(:);

n_nodes    = length(x);
n_internal = n_nodes - 2;

%% ========================================================================
% Allocate
%% ========================================================================

a_WW = zeros(n_internal,1);
a_W  = zeros(n_internal,1);
a_P  = zeros(n_internal,1);
a_E  = zeros(n_internal,1);
a_EE = zeros(n_internal,1);

rhs = zeros(n_internal,1);

has_pentadiagonal = false;

%% ========================================================================
% Assembly
%% ========================================================================

for I = 1:n_internal

    i = I + 1;

    is_new_format = false;

    try

        [coeffs,offsets,source] = ...
            scheme(mesh,U,k,i);

        if isnumeric(offsets) && ...
           all(ismember(offsets(:)',[-2 -1 0 1 2]))

            is_new_format = true;

        end

    catch

        is_new_format = false;

    end

    if is_new_format

        for m = 1:length(offsets)

            switch offsets(m)

                case -2
                    a_WW(I)=coeffs(m);

                case -1
                    a_W(I)=coeffs(m);

                case 0
                    a_P(I)=coeffs(m);

                case 1
                    a_E(I)=coeffs(m);

                case 2
                    a_EE(I)=coeffs(m);

                otherwise
                    error('Unsupported stencil offset.');

            end

        end

        if any(abs(offsets)==2)

            has_pentadiagonal = true;

        end

    else

        [aW,aP,aE,source] = ...
            scheme(mesh,U,k,i);

        a_W(I)=aW;
        a_P(I)=aP;
        a_E(I)=aE;

    end

    rhs(I)=source;

end

%% ========================================================================
% Periodic BCs
%% ========================================================================

is_periodic = ...
    strcmpi(bc.W.type,'periodic') && ...
    strcmpi(bc.E.type,'periodic');

if is_periodic

    %
    % Do nothing.
    %
    % Periodic couplings are handled later
    % in the transient solver / matrix assembly.
    %

%% ========================================================================
% Dirichlet / Neumann BCs
%% ========================================================================

elseif strcmpi(bc.W.type,'dirichlet')

    rhs(1) = rhs(1) ...
           - a_W(1)*bc.W.value;

    if has_pentadiagonal && n_internal >= 2

        rhs(2) = rhs(2) ...
               - a_WW(2)*bc.W.value;

    end

elseif strcmpi(bc.W.type,'neumann')

    dxW = x(2)-x(1);

    a_P(1) = a_P(1)+a_W(1);

    rhs(1) = rhs(1) ...
           + a_W(1)*bc.W.value*dxW;

else

    error('Unknown west BC.');

end

%% East BC

if ~is_periodic

    if strcmpi(bc.E.type,'dirichlet')

        rhs(end)=rhs(end) ...
               - a_E(end)*bc.E.value;

        if has_pentadiagonal && n_internal >= 2

            rhs(end-1)=rhs(end-1) ...
                      - a_EE(end-1)*bc.E.value;

        end

    elseif strcmpi(bc.E.type,'neumann')

        dxE = x(end)-x(end-1);

        a_P(end)=a_P(end)+a_E(end);

        rhs(end)=rhs(end) ...
                - a_E(end)*bc.E.value*dxE;

    else

        error('Unknown east BC.');

    end

end

%% ========================================================================
% Store
%% ========================================================================

A.a_WW = a_WW;
A.a_W  = a_W;
A.a_P  = a_P;
A.a_E  = a_E;
A.a_EE = a_EE;

A.periodic = is_periodic;

end