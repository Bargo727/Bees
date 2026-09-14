% Figure 1: Baseline Time Series
clear; clc; close all;

% Parameters
p.alpha = 1.5e-9;   % Bee growth rate scalar
p.K = 50000;      % Carrying capacity
p.A = 1000;       % Allee threshold
p.beta = 5e-5;    % Bee death rate due to phoretic mites
p.nu = 1.5;       % Mite reproduction factor
p.sigma = 0.077;    % Brood emergence rate
p.gamma = 2e-6;   % Brood invasion rate
p.delta = 0.035;   % Natural death of phoretic mites

tspan = [0 6000]; % 150 days

B_star = p.delta/(p.gamma*(p.nu-1));

% Initial Conditions: [Bees, Phoretic Mites, Repro Mites]
ICs = [
    8000, 100,  50;    % Case 1: Endemic Trap
    900,  10,   0;     % Case 2: Allee Collapse
    9000, 4000, 2000   % Case 3: No parasite
];

colors = {'#77AC30', '#D95319', '#000000'};
labels = {'endemic trap', 'Allee collapse', 'parasite overwhelm'};

figure('Position', [100, 100, 800, 400]);
subplot(1,2,1); hold on; grid on;
subplot(1,2,2); hold on; grid on;

for i = 1:size(ICs, 1)
    [t, y] = ode45(@(t,y) bee_mite_ode(t, y, p), tspan, ICs(i,:));
    
    subplot(1,2,1); % Bee Population
    plot(t, y(:,1), 'Color', colors{i}, 'LineWidth', 2, 'DisplayName', labels{i});
    
    subplot(1,2,2); % Total Mite Population (Phoretic + Repro)
    plot(t, y(:,2) + y(:,3), 'Color', colors{i}, 'LineWidth', 2, 'DisplayName', labels{i});
    box off
end

subplot(1,2,1);
yline(p.A, 'r--', 'Allee threshold', 'HandleVisibility', 'off');
yline(p.K, 'k--', 'carrying capacity', 'HandleVisibility', 'off');
xlabel('time (days)'); ylabel('bee population');
title('Bee Dynamics'); legend('Location', 'best');
set(gca,'fontsize',20)

subplot(1,2,2);
xlabel('time (days)'); ylabel('total mite population');
title('Mite Dynamics');
set(gca,'fontsize',20)

function dydt = bee_mite_ode(~, y, p)
    B = y(1); Mp = y(2); Mr = y(3);
    dBdt = p.alpha * B * (B - p.A) * (p.K - B) - p.beta * B * Mp;
    dMpdt = p.nu * p.sigma * Mr - p.gamma * B * Mp - p.delta * Mp;
    dMrdt = -p.sigma * Mr + p.gamma * B * Mp;
    dydt = [dBdt; dMpdt; dMrdt];
end