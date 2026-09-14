% Figure: Multi-Panel Dynamical Regimes of the Bee-Mite Model
clear; clc; close all;

% --- Base Parameters ---
p.alpha = 1.5e-9;   % Bee growth rate scalar
p.A = 1000;         % Allee threshold
p.beta = 5e-5;      % Bee death rate due to phoretic mites
p.nu = 1.5;         % Mite reproduction factor
p.sigma = 0.077;    % Brood emergence rate
p.gamma = 2e-6;     % Brood invasion rate
p.delta = 0.035;    % Natural death of phoretic mites

% Calculate B^* (Critical bee population for mites)
B_star = p.delta / (p.gamma * (p.nu - 1)); % Equals 35,000

tspan = [0 800]; % 800 days is enough to show asymptotic convergence
opts = odeset('NonNegative', 1:3); % Prevents numerical overshoot below zero

% --- Define the Three Scenarios ---
% Scenario 1: Zero Mites
S(1).K = 50000;
S(1).ICs = [
    8000, 0, 0;  % Recovers to K
    900,  0, 0   % Allee Collapse
];
S(1).colors = {'#77AC30', '#D95319'};
S(1).labels = {'healthy hive', 'extinction'};
S(1).title_bee = 'No Mites';
S(1).title_mite = 'No Mites';

% Scenario 2: K < B^* (Healthy Hive)
S(2).K = 30000; % K is below B^* = 35000
S(2).ICs = [
    8000, 1000, 500;  % Mites die out, Bees recover to K
    900,  100,  50    % Allee Collapse
];
S(2).colors = {'#77AC30', '#D95319'};
S(2).labels = {'bee recovery', 'extinction'};
S(2).title_bee = 'K < B^* (Healthy Hive)';
S(2).title_mite = 'K < B^* (Mites Die Out)';

% Scenario 3: K > B^* (Endemic / Overwhelm)
S(3).K = 50000; % K is above B^* = 35000
S(3).ICs = [
    38000, 10000, 5000;   % Endemic Trap
    900,   10,    0;      % Allee Collapse
    45000, 35000, 10000   % Parasite Overwhelm
];
S(3).colors = {'#0072BD', '#D95319', '#000000'};
S(3).labels = {'Endemic Trap', 'Allee Collapse', 'Parasite Overwhelm'};
S(3).title_bee = 'K > B^* (Endemicity & Collapse)';
S(3).title_mite = 'K > B^* (Mite Dynamics)';

% --- Plotting ---
figure('Position', [100, 100, 1000, 1200]);


for row = 1:3
    p.K = S(row).K; % Update carrying capacity for this regime
    
    % Bee Subplot
    ax_bee = subplot(3, 2, 2*row - 1); hold on; grid on;
    % Mite Subplot
    ax_mite = subplot(3, 2, 2*row); hold on; grid on;
    
    for i = 1:size(S(row).ICs, 1)
        [t, y] = ode45(@(t,y) bee_mite_ode(t, y, p), tspan, S(row).ICs(i,:), opts);
        
        % Plot Bees
        plot(ax_bee, t, y(:,1), 'Color', S(row).colors{i}, 'LineWidth', 4, 'DisplayName', S(row).labels{i});
        
        % Plot Total Mites (Phoretic + Repro)
        plot(ax_mite, t, y(:,2) + y(:,3), 'Color', S(row).colors{i}, 'LineWidth', 4, 'DisplayName', S(row).labels{i});
    end
    
    % Format Bee Plot
    yline(ax_bee, p.K, 'k--', 'K', 'HandleVisibility', 'off', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    yline(ax_bee, p.A, 'r--', 'A', 'HandleVisibility', 'off', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    yline(ax_bee, B_star, 'b--', 'B^*', 'HandleVisibility', 'off', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    
    xlabel(ax_bee, 'Time (days)'); 
    ylabel(ax_bee, 'Bee Population');
    title(ax_bee, S(row).title_bee); 
    legend(ax_bee, 'Location', 'best');
    set(ax_bee, 'fontsize', 12, 'YLim', [0 55000]);
    set(gca,'fontsize',20)
    
    % Format Mite Plot
    xlabel(ax_mite, 'Time (days)'); 
    ylabel(ax_mite, 'Total Mite Population');
    title(ax_mite, S(row).title_mite);
    set(ax_mite, 'fontsize', 12);
    set(gca,'fontsize',20)
    if row == 1
        set(ax_mite, 'YLim', [0 100]); % Keep axis scaled nicely when zero
    end
end

% --- ODE Function ---
function dydt = bee_mite_ode(~, y, p)
    B = y(1); Mp = y(2); Mr = y(3);
    
    % Force derivatives to 0 if bees collapse to prevent numerical artifacts
    if B <= 1e-3
        dBdt = 0;
        dMpdt = -p.delta * Mp;
        dMrdt = -p.sigma * Mr;
    else
        dBdt = p.alpha * B * (B - p.A) * (p.K - B) - p.beta * B * Mp;
        dMpdt = p.nu * p.sigma * Mr - p.gamma * B * Mp - p.delta * Mp;
        dMrdt = -p.sigma * Mr + p.gamma * B * Mp;
    end
    
    dydt = [dBdt; dMpdt; dMrdt];
end