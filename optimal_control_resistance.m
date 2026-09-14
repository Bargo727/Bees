function optimal_control_resistance()
    % =========================================================
    % 1. Parameter Definitions (Scaled to 10k as in previous models)
    % =========================================================
    T = 350;                % Final time (days)
    N = 2000;               % Number of time steps
    dt = T/N;               % Time step size
    t = linspace(0, T, N+1);
    
    % Biological Parameters
    alpha = 2E-8; A = 1000; K = 10000.0; 
    beta = 2E-5;             % true_beta (5e-5) * 10000
    eta = 0.008; nu = 2.0; sigma = 0.1; 
    delta = 0.05; gamma = 1E-5;
    
    % Resistance Parameters
    mu = 0.02;              % Mutation probability
    rho = 0.1;              % Efficacy of u2 on resistant mites (10%)
    
    % Objective Functional Weights (Requires tuning)
    w1 = 1000.0;             % Payoff for final bee population
    w2 = 0.1;               % Penalty for phoretic mites
    w3 = 0.1;               % Penalty for reproductive mites
    C1 = 10.0;              % Cost weight for soft treatment (u1)
    C2 = 400.0;              % Cost weight for harsh treatment (u2)
    
    u1_max = 1.0; u2_max = 1.0; % Upper bounds on controls
    
    % =========================================================
    % 2. Initialization
    % =========================================================
    x = zeros(5, N+1);      % States: [B; Mps; Mrs; Mpr; Mrr]
    l = zeros(5, N+1);      % Adjoints: [L1; L2; L3; L4; L5]
    u1 = zeros(1, N+1);     % Control 1 (Soft)
    u2 = zeros(1, N+1);     % Control 2 (Harsh)
    
    % Initial Conditions
    x(:, 1) = [8000; 5000; 2500; 1000; 500]; 
    
    tol = 1e-4;             % Convergence tolerance
    max_iter = 1000;        % Maximum FBSM iterations
    test = -1; iter = 0;
    
    % =========================================================
    % 3. Forward-Backward Sweep Loop
    % =========================================================
    while (test < 0) && (iter < max_iter)
        iter = iter + 1;
        
        old_u1 = u1; old_u2 = u2;
        old_x = x; old_l = l;
        
        % --- FORWARD SWEEP (State Equations via RK4) ---
        for i = 1:N
            u1_mid = 0.5 * (u1(i) + u1(i+1));
            u2_mid = 0.5 * (u2(i) + u2(i+1));
            
            k1 = state_eqs(x(:, i), u1(i), u2(i));
            k2 = state_eqs(x(:, i) + 0.5*dt*k1, u1_mid, u2_mid);
            k3 = state_eqs(x(:, i) + 0.5*dt*k2, u1_mid, u2_mid);
            k4 = state_eqs(x(:, i) + dt*k3, u1(i+1), u2(i+1));
            
            x(:, i+1) = x(:, i) + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        end
        
        % --- BACKWARD SWEEP (Adjoint Equations via RK4) ---
        % Transversality conditions
        l(:, N+1) = [w1; 0; 0; 0; 0]; 
        
        for i = N:-1:1
            u1_mid = 0.5 * (u1(i) + u1(i+1));
            u2_mid = 0.5 * (u2(i) + u2(i+1));
            x_mid = 0.5 * (x(:, i) + x(:, i+1));
            
            k1 = adjoint_eqs(l(:, i+1), x(:, i+1), u1(i+1), u2(i+1));
            k2 = adjoint_eqs(l(:, i+1) - 0.5*dt*k1, x_mid, u1_mid, u2_mid);
            k3 = adjoint_eqs(l(:, i+1) - 0.5*dt*k2, x_mid, u1_mid, u2_mid);
            k4 = adjoint_eqs(l(:, i+1) - dt*k3, x(:, i), u1(i), u2(i));
            
            l(:, i) = l(:, i+1) - (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
        end
        
        % --- UPDATE CONTROLS ---
        % Calculate optimal characterization
        u1_new = zeros(1, N+1);
        u2_new = zeros(1, N+1);
        
        for i = 1:N+1
            B = x(1,i); Mps = x(2,i); Mrs = x(3,i); Mrr = x(5,i);
            L1 = l(1,i); L2 = l(2,i); L3 = l(3,i); L5 = l(5,i);
            
            % Pontryagin characterizations
            val1 = (-L2 * Mps) / (2 * C1);
            val2 = (-L1 * eta * B - L3 * Mrs - L5 * rho * Mrr) / (2 * C2);
            
            u1_new(i) = min(u1_max, max(0, val1));
            u2_new(i) = min(u2_max, max(0, val2));
        end
        
        % Convex combination to ensure stability
        u1 = 0.5 * u1_new + 0.5 * u1;
        u2 = 0.5 * u2_new + 0.5 * u2;
        
        % Check Convergence
        err1 = sum(abs(u1 - old_u1)) / sum(abs(u1) + eps);
        err2 = sum(abs(u2 - old_u2)) / sum(abs(u2) + eps);
        err3 = sum(sum(abs(x - old_x))) / sum(sum(abs(x) + eps));
        
        if max([err1, err2, err3]) < tol
            test = 1;
            fprintf('Converged in %d iterations.\n', iter);
        end
    end

    % =========================================================
    % 4. Plotting Results
    % =========================================================
    figure;
    subplot(2,1,1); hold on; grid on;
    plot(t, u1, 'b', 'LineWidth', 4); hold on;
    plot(t, u2, 'r', 'LineWidth', 4);
    title('Optimal Treatment Cocktail with Resistant Mutants');
    h = legend('soft treatment', 'harsh treatment');
    set(h,'box','off')
    ylabel('treatment intensity');
    set(gca,'fontsize',20)
    xlim([0 340])
    
    subplot(2,1,2); hold on; grid on;
    plot(t, x(1,:), 'g', 'LineWidth', 4); hold on;
    plot(t, sum(x(2:5,:)), 'k--', 'LineWidth', 4);
    %title('Population Dynamics');
    h = legend('bee population', 'total mites');
    set(h,'box','off')
    ylabel('populations');
    xlabel('time (days)');
    set(gca,'fontsize',20)
    xlim([0 340])

    % =========================================================
    % Nested Functions for System Dynamics
    % =========================================================
    function dx = state_eqs(xx, u1_val, u2_val)
        B = xx(1); Mps = xx(2); Mrs = xx(3); Mpr = xx(4); Mrr = xx(5);
        dx = zeros(5,1);
        
        dx(1) = alpha*B*(B-A)*(K-B) - beta*B*(Mps + Mpr) - eta*u2_val*B;
        dx(2) = (1-mu)*nu*sigma*Mrs - gamma*B*Mps - (delta + u1_val)*Mps;
        dx(3) = -(sigma + u2_val)*Mrs + gamma*B*Mps;
        dx(4) = nu*sigma*Mrr + mu*nu*sigma*Mrs - gamma*B*Mpr - delta*Mpr;
        dx(5) = -(sigma + rho*u2_val)*Mrr + gamma*B*Mpr;
    end

    function dl = adjoint_eqs(ll, xx, u1_val, u2_val)
        B = xx(1); Mps = xx(2); Mrs = xx(3); Mpr = xx(4); Mrr = xx(5);
        l1 = ll(1); l2 = ll(2); l3 = ll(3); l4 = ll(4); l5 = ll(5);
        dl = zeros(5,1);
        
        dl(1) = -l1*(alpha*(-3*B^2 + 2*(A+K)*B - A*K) - beta*(Mps + Mpr) - eta*u2_val) ...
                + gamma*Mps*(l2 - l3) + gamma*Mpr*(l4 - l5);
        dl(2) = w2 + l1*beta*B + l2*(gamma*B + delta + u1_val) - l3*gamma*B;
        dl(3) = w3 - l2*(1-mu)*nu*sigma - l4*mu*nu*sigma + l3*(sigma + u2_val);
        dl(4) = w2 + l1*beta*B + l4*(gamma*B + delta) - l5*gamma*B;
        dl(5) = w3 - l4*nu*sigma + l5*(sigma + rho*u2_val);
    end
end