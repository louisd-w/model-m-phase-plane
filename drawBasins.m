
function drawBasins(ax, p, xmax, ymax)

    % unpack parameters
    a   = p.a;
    b   = p.b;
    c   = p.c;
    d   = p.d;
    phi = p.phi;
    K   = p.K;
    a_f = p.a_f;
    W_c = p.W_c;
    K_b = p.K_b;
    basin_grid_size = p.basin_grid_size;
    
    % define functions
    mu = @(M) d + phi*c.*M;
    beta = @(M) (a.*M)./(b+M) - (1-phi)*c.*M;
    U = @(M) a_f.*M + (W_c.*M)./(K_b+M);
    U_prime = @(M) a_f + (W_c*K_b)./(K_b+M).^2;
    S = @(C,M) C .* (beta(M).*(1-C./K) - mu(M));
    
    % Model M
    F = @(C,M) S(C,M);
    G = @(C,M) -U(M).*S(C,M)./(1 + C.*U_prime(M));
    f = @(t,x) [F(x(1),x(2)); G(x(1),x(2))];
    
    % time interval for each trajectory
    T = [0, 100];
    
    % make grid for basin plot, sampling from the centres
    dC = xmax/basin_grid_size;
    dM = ymax/basin_grid_size;
    C_vals = linspace(dC/2, xmax-dC/2, basin_grid_size);
    M_vals = linspace(dM/2, ymax-dM/2, basin_grid_size);
    [C_grid,M_grid] = meshgrid(C_vals,M_vals);
   
    % store classification:
        % 0 = unclassified
        % 1 = extinction
        % 2 = positive equilibrium
    basin = zeros(size(C_grid));

    % solve from each grid point, x = (C,M)
    for i = 1:numel(C_grid)

        x0 = [C_grid(i); M_grid(i)];

        [t_sol, x_sol] = ode45(f, T, x0);

        % final point of trajectory
        Cf = x_sol(end,1);
        Mf = x_sol(end,2);
        
        % distance from extinction equilibrium C = 0
        d_ext = abs(Cf);
        
        % check whether a positive equilibrium exists at Mf
        if abs(beta(Mf)) < 1e-10
        
            % beta is effectively zero, so positive equilibrium formula is invalid
            basin(i) = 1;
        
        else
        
            % positive equilibrium at M = Mf
            Ceq = K * (1 - mu(Mf)/beta(Mf));
        
            if Ceq <= 0
        
                % no positive equilibrium at this value of M
                basin(i) = 1;
        
            else
        
                % distance from positive equilibrium
                d_pos = abs(Cf - Ceq);
        
                % decide which equilibrium branch is closer
                if d_ext < d_pos
                    basin(i) = 1;
                else
                    basin(i) = 2;
                end
        
            end
        
        end
    end
    
    % clear plot
    cla(ax)
    hold(ax, 'on')

    % labeling and sizing plot
    xlabel(ax, 'C')
    ylabel(ax, 'M')
    xlim(ax, [0 xmax])
    ylim(ax, [0 ymax])
    daspect(ax, [1 1 1]) % keep axis in proportion
    grid(ax, 'on')
    title(ax, 'Click initial position for a trajectory')

    % shade basins
    h = imagesc(ax, C_grid(1,:), M_grid(:,1), basin);
    
    % basin colours
    colormap(ax, [
        0 0 1          % extinction
        1 0 0     % positive equilibrium
    ])
    
    clim(ax, [1 2])
    
    % transparent where unclassified and partially transparent elsewhere
    h.AlphaData = 0.25*(basin ~= 0);
    
    % turn off mouse
    h.HitTest = 'off';
    h.PickableParts = 'none';
    
    hold(ax, 'off')
    
end