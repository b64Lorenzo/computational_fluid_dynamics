function [phi,time,A,rhs] = ...
    unsteady_convection_diffusion_1d( ...
    mesh,U,k,bc,scheme,...
    dt,t_end,time_scheme,fid)

if nargin < 9
    fid = [];
end

if nargin < 8 || isempty(time_scheme)
    time_scheme = 'implicit_euler';
end

%% ========================================================================
% Mesh
%% ========================================================================

x = mesh.x(:);

n_nodes    = length(x);
n_internal = n_nodes - 2;

time = 0:dt:t_end;
n_steps = length(time);

%% ========================================================================
% Initial condition
%% ========================================================================

phi = zeros(n_nodes,n_steps);

if isfield(bc,'IC')

    if length(bc.IC) ~= n_nodes

        error( ...
            'IC length (%d) differs from mesh size (%d).', ...
            length(bc.IC), ...
            n_nodes);

    end

    phi(:,1) = bc.IC(:);

end

%% ========================================================================
% Spatial operator
%% ========================================================================

[A,rhs,has_pentadiagonal] = ...
    assemble_convection_diffusion_1d( ...
    mesh,U,k,bc,scheme);

a_WW = A.a_WW;
a_W  = A.a_W;
a_P  = A.a_P;
a_E  = A.a_E;
a_EE = A.a_EE;

is_periodic = false;

if isfield(A,'periodic')
    is_periodic = A.periodic;
end

%% ========================================================================
% Time loop
%% ========================================================================

for n = 1:n_steps-1

    phi_old = phi(2:end-1,n);

    switch lower(time_scheme)

        %==============================================================
        % Explicit Euler
        %==============================================================

        case {'explicit','explicit_euler'}

            residual = zeros(n_internal,1);

            for I = 1:n_internal

                if is_periodic

                    iW  = mod(I-2,n_internal)+1;
                    iE  = mod(I,n_internal)+1;

                    r = ...
                        a_W(I)*phi_old(iW) + ...
                        a_P(I)*phi_old(I)  + ...
                        a_E(I)*phi_old(iE);

                    if has_pentadiagonal

                        iWW = mod(I-3,n_internal)+1;
                        iEE = mod(I+1,n_internal)+1;

                        r = r ...
                          + a_WW(I)*phi_old(iWW) ...
                          + a_EE(I)*phi_old(iEE);

                    end

                else

                    r = a_P(I)*phi_old(I);

                    if I > 1
                        r = r + a_W(I)*phi_old(I-1);
                    end

                    if I < n_internal
                        r = r + a_E(I)*phi_old(I+1);
                    end

                    if has_pentadiagonal

                        if I > 2
                            r = r + a_WW(I)*phi_old(I-2);
                        end

                        if I < n_internal-1
                            r = r + a_EE(I)*phi_old(I+2);
                        end

                    end

                end

                residual(I) = r - rhs(I);

            end

            %
            % Validated sign convention
            %

            phi_new = phi_old + dt*residual;

        %==============================================================
        % Implicit Euler
        %==============================================================

        case {'implicit','implicit_euler'}

            rhs_t = phi_old/dt - rhs;

            if is_periodic

                M = sparse(n_internal,n_internal);

                for I = 1:n_internal

                    M(I,I) = 1/dt - a_P(I);

                    iW = mod(I-2,n_internal)+1;
                    iE = mod(I,n_internal)+1;

                    M(I,iW) = M(I,iW) - a_W(I);
                    M(I,iE) = M(I,iE) - a_E(I);

                    if has_pentadiagonal

                        iWW = mod(I-3,n_internal)+1;
                        iEE = mod(I+1,n_internal)+1;

                        M(I,iWW) = ...
                            M(I,iWW) - a_WW(I);

                        M(I,iEE) = ...
                            M(I,iEE) - a_EE(I);

                    end

                end

                phi_new = M\rhs_t;

            else

                aP_t = 1/dt - a_P;

                if has_pentadiagonal

                    b2 = zeros(n_internal,1);
                    b1 = zeros(n_internal,1);
                    c1 = zeros(n_internal,1);
                    c2 = zeros(n_internal,1);

                    b2(3:end)   = -a_WW(3:end);
                    b1(2:end)   = -a_W(2:end);

                    c1(1:end-1) = -a_E(1:end-1);
                    c2(1:end-2) = -a_EE(1:end-2);

                    phi_new = pdma( ...
                        b2,...
                        b1,...
                        aP_t,...
                        c1,...
                        c2,...
                        rhs_t,...
                        n_internal);

                else

                    lower_diag = zeros(n_internal,1);
                    upper_diag = zeros(n_internal,1);

                    lower_diag(2:end)   = -a_W(2:end);
                    upper_diag(1:end-1) = -a_E(1:end-1);

                    phi_new = tdma( ...
                        lower_diag,...
                        aP_t,...
                        upper_diag,...
                        rhs_t,...
                        n_internal);

                end

            end

        otherwise

            error('Unknown time scheme.');

    end

    %% ====================================================================
    % Store interior
    %% ====================================================================

    phi(2:end-1,n+1) = phi_new;

    %% ====================================================================
    % Boundary treatment
    %% ====================================================================

    if is_periodic

        phi(1,n+1)   = phi(end-1,n+1);
        phi(end,n+1) = phi(2,n+1);

    else

        %--------------------------------------------------------------
        % West
        %--------------------------------------------------------------

        if strcmpi(bc.W.type,'dirichlet')

            phi(1,n+1) = bc.W.value;

        elseif strcmpi(bc.W.type,'neumann')

            dxW = x(2)-x(1);

            phi(1,n+1) = ...
                phi(2,n+1) - bc.W.value*dxW;

        end

        %--------------------------------------------------------------
        % East
        %--------------------------------------------------------------

        if strcmpi(bc.E.type,'dirichlet')

            phi(end,n+1) = bc.E.value;

        elseif strcmpi(bc.E.type,'neumann')

            dxE = x(end)-x(end-1);

            phi(end,n+1) = ...
                phi(end-1,n+1) + bc.E.value*dxE;

        end

    end

    %% ====================================================================
    % Safety check
    %% ====================================================================

    if any(~isfinite(phi_new))

        error( ...
            'NaN/Inf detected at timestep %d', ...
            n);

    end

end

end