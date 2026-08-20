
function drawTrajectory(ax, p, C0, M0)

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

    % time interval for each trajectory
    T = linspace(0,100,5000);

    % solve from clicked point, x = (C,M)
    x0 = [C0; M0];
    f = @(t,x) [F(x(1),x(2)); G(x(1),x(2))];
    [t_sol, x_sol] = ode45(f, T, x0);

    % plot initial point
    plot(ax, C0, M0, 'o', 'MarkerFaceColor', 'g')

    % plot trajectory
    h = plot(ax, x_sol(:,1), x_sol(:,2), 'g-', 'LineWidth', 1.5);
    h.HitTest = 'off';
    h.PickableParts = 'none';

end