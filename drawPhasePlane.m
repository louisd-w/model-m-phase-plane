
function drawPhasePlane(ax, p, xmax, ymax)

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

    % define functions
    mu = @(M) d + phi*c.*M;
    beta = @(M) (a.*M)./(b+M) - (1-phi)*c.*M;
    U = @(M) a_f.*M + (W_c.*M)./(K_b+M);
    U_prime = @(M) a_f + (W_c*K_b)./(K_b+M).^2;
    S = @(C,M) C .* (beta(M).*(1-C./K) - mu(M));

    % Model M
    F = @(C,M) S(C,M);
    G = @(C,M) -U(M).*S(C,M)./(1 + C.*U_prime(M));

    % grid for vector field
    [C_grid,M_grid] = meshgrid(linspace(0,xmax,20), linspace(0,ymax,20));

    % evaluate vector field
    dC = F(C_grid,M_grid);
    dM = G(C_grid,M_grid);

    % normalise vectors
    L = sqrt(dC.^2 + dM.^2);
    L(L == 0) = 1;
    dC = dC ./ L;
    dM = dM ./ L;

    % clear plot
    cla(ax)
    hold(ax, 'on')

    % vector field
    q = quiver(ax, C_grid, M_grid, dC, dM, 0.5);
    % turn off mouse
    q.HitTest = 'off';
    q.PickableParts = 'none';

    % nullclines
    f1 = fimplicit(ax, F, [0 xmax 0 ymax], '--', 'Color', [0 0 1], 'LineWidth', 2.5);
    f2 = fimplicit(ax, G, [0 xmax 0 ymax], 'r-', 'LineWidth', 1.5);
    % turn off mouse
    f1.HitTest = 'off';
    f1.PickableParts = 'none';
    f2.HitTest = 'off';
    f2.PickableParts = 'none';
    % legend
    lgd = legend(ax, [f1 f2], {'$\dot{C}=0$', '$\dot{M}=0$'}, 'Interpreter', 'latex', 'Location', 'northeast');
    lgd.AutoUpdate = 'off';

    % labeling and sizing plot
    xlabel(ax, 'C')
    ylabel(ax, 'M')
    xlim(ax, [0 xmax])
    ylim(ax, [0 ymax])
    daspect(ax, [1 1 1]) % keep axis in proportion
    grid(ax, 'on')
    title(ax, 'Click initial position for a trajectory')

    % make phase plane respond reliably to clicks
    ax.HitTest = 'on';
    ax.PickableParts = 'all';
    disableDefaultInteractivity(ax)
    ax.Toolbar = [];

end











