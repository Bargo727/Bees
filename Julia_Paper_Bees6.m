% Bulletproof Optimal Control of the Treatment Cocktail
clear; clc; close all;

% --- 1. Parameters ---
p.alpha = 2e-8; p.K = 10000; p.A = 1000; 
p.beta  = 2e-5; % Lowered lethality so it's a chronic drain, not a sudden execution
p.nu = 2.0; p.sigma = 0.1; 
p.gamma = 1e-5; p.delta = 0.05; 
p.eta = 0.008;  % Lowered toxicity so the solver isn't terrified to use u2

% Objective Function Weights (The "Rescue" Economics)
w1 = 1000;   % Massive reward for final bees (SAVE THEM AT ALL COSTS)
w2 = 0.1;    % Penalty for phoretic mites
w3 = 0.1;    % Penalty for reproductive mites
C1 = 10;      % Dirt cheap Soft Treatment
C2 = 400;     % Dirt cheap Harsh Treatment

% Simulation Setup
T = 350;            
dt = 0.1;           
time = 0:dt:T;
N = length(time);

% Initial Conditions
B0 = 8000; Mp0 = 5000; Mr0 = 2500;

% --- 2. Baseline Simulation (No Control) ---
[t_base, y_base] = ode45(@(t,y) [ ...
    p.alpha*y(1)*(y(1)-p.A)*(p.K-y(1)) - p.beta*y(1)*y(2); ...
    p.nu*p.sigma*y(3) - p.gamma*y(1)*y(2) - p.delta*y(2); ...
    -p.sigma*y(3) + p.gamma*y(1)*y(2)], time, [B0, Mp0, Mr0]);

% --- 3. Optimal Control (Forward-Backward Sweep) ---
B = zeros(1, N); Mp = zeros(1, N); Mr = zeros(1, N);
B(1) = B0; Mp(1) = Mp0; Mr(1) = Mr0;

L1 = zeros(1, N); L2 = zeros(1, N); L3 = zeros(1, N);
u1 = ones(1, N); u2 = ones(1, N);

max_iter = 2000; 
tolerance = 1e-3; 
test = 1; 
iter = 0;

% *** CRITICAL FIX: The Relaxation Parameter ***
update_weight = 0.05; % Only take a 5% step toward the new optimal control

while (test > tolerance && iter < max_iter)
    iter = iter + 1;
    
    oldu1 = u1; oldu2 = u2;
    
    % FORWARD SWEEP
    for i = 1:N-1
        dB = p.alpha*B(i)*(B(i)-p.A)*(p.K-B(i)) - p.beta*B(i)*Mp(i) - p.eta*u2(i)*B(i);
        dMp = p.nu*p.sigma*Mr(i) - p.gamma*B(i)*Mp(i) - (p.delta + u1(i))*Mp(i);
        dMr = -(p.sigma + u2(i))*Mr(i) + p.gamma*B(i)*Mp(i);
        
        B(i+1)  = max(0, B(i)  + dt * dB);
        Mp(i+1) = max(0, Mp(i) + dt * dMp);
        Mr(i+1) = max(0, Mr(i) + dt * dMr);
    end
    
    % BACKWARD SWEEP
    L1(N) = w1; L2(N) = 0; L3(N) = 0;
    
    for i = N:-1:2
        dH_dB = L1(i)*(p.alpha*(-3*B(i)^2 + 2*(p.A+p.K)*B(i) - p.A*p.K) - p.beta*Mp(i) - p.eta*u2(i)) ...
              + L2(i)*(-p.gamma*Mp(i)) + L3(i)*(p.gamma*Mp(i));
          
        dH_dMp = -w2 + L1(i)*(-p.beta*B(i)) + L2(i)*(-p.gamma*B(i) - p.delta - u1(i)) + L3(i)*(p.gamma*B(i));
        dH_dMr = -w3 + L2(i)*(p.nu*p.sigma) + L3(i)*(-p.sigma - u2(i));
        
        L1(i-1) = L1(i) + dt * dH_dB;
        L2(i-1) = L2(i) + dt * dH_dMp;
        L3(i-1) = L3(i) + dt * dH_dMr;
    end
    
    % UPDATE CONTROLS
    for i = 1:N
        temp_u1 = (-L2(i) * Mp(i)) / (2 * C1);
        temp_u2 = (-L3(i) * Mr(i) - L1(i) * p.eta * B(i)) / (2 * C2);
        
        temp_u1 = max(0, min(0.6, temp_u1)); % Max soft dose
        temp_u2 = max(0, min(0.4, temp_u2)); % Max harsh dose
        
        % Slow, stable update to prevent falling off the Allee cliff
        u1(i) = update_weight * temp_u1 + (1 - update_weight) * oldu1(i);
        u2(i) = update_weight * temp_u2 + (1 - update_weight) * oldu2(i);
    end
    
    test = sum(abs(u1 - oldu1)) + sum(abs(u2 - oldu2));
end

disp(['Converged in ', num2str(iter), ' iterations.']);

% --- 4. Plotting ---
figure('Position', [300, 50, 700, 800]);

subplot(3,1,1); hold on; grid on;
plot(time, u1, 'b-', 'LineWidth', 4, 'DisplayName', 'optimal soft treatment');
plot(time, u2, 'r-', 'LineWidth', 4, 'DisplayName', 'optimal harsh treatment');
ylabel('treatment intensity');
title('Optimal Control Strategy', 'FontSize', 14);
h = legend('Location', 'northeast');
set(h,'box','off')
xlim([0 340])
set(gca,'fontsize',20)

subplot(3,1,2); hold on; grid on;
plot(time, Mp, 'Color', '#D95319', 'LineWidth', 3, 'DisplayName', 'phoretic (optimal)');
plot(time, Mr, 'Color', '#EDB120', 'LineWidth', 3, 'DisplayName', 'reproductive (optimal)');
plot(t_base, y_base(:,2)+y_base(:,3), 'k--', 'LineWidth', 1.5, 'DisplayName', 'total mites (untreated)');
ylabel('mite population');
%title('Parasite Suppression', 'FontSize', 14);
h = legend('Location', 'northeast');
set(h,'box','off')
xlim([0 340])
set(gca,'fontsize',20)

subplot(3,1,3); hold on; grid on;
plot(time, B, 'g-', 'LineWidth', 3, 'DisplayName', 'optimal cocktail survival');
plot(t_base, y_base(:,1), 'k--', 'LineWidth', 2, 'DisplayName', 'untreated collapse');
yline(p.A, 'r:', 'LineWidth', 2, 'DisplayName', 'Allee threshold');
xlabel('time (days)');
ylabel('bee population');
%title('Colony Rescue', 'FontSize', 14);
h = legend('Location', 'east');
set(h,'box','off')
set(gca,'fontsize',20)
xlim([0 340])
ylim([0 11000]);