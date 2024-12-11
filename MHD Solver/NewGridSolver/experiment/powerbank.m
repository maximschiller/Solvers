clc
clear

clear all
close all

% power bank parameters
mu = (10^-7)*4*pi();
ra = 10/1000;
rc = 5/1000;
ncapsseries = 5;
ncapsparallel = 5;
numcaps = ncapsseries * ncapsparallel;
V = 1350;
Vcap = 400;
IndivCcap = 10e-6;

Ccap = ncapsparallel * (1 / (ncapsseries / (IndivCcap)));
Rres = 1 / (ncapsparallel / (ncapsseries * 0.1));
Lind = 1 / (ncapsparallel / (ncapsseries * (10e-9)));

disp(['Capacitance: ' num2str(Ccap)]);
disp(['Resistance = ' num2str(Rres)]);
disp(['Inductance = ' num2str(Lind)]);

RLsq = (Rres / (2 * Lind))^2;
LC = 1 / (Lind * Ccap);
disp(['(R/L)^2 = ' num2str(RLsq) ' LC = ' num2str(LC)]); % need to ensure RLsq is greater than LC

% Plot the discharge current vs time
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9); % optional
[t, sol] = ode45(@f, [0 10e-6], [0 0], options);

charge = sol(:, 1);
current = sol(:, 2);

% Plot the discharge current vs time
figure;
plot(t, current, 'blue', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Current (A)');
% set(gca, 'TickLabelInterpreter', 'latex');
% set(gca, 'TickLabelFormat', '%.1e');
grid on;
grid minor
box on
% saveas(gcf, 'dischargeA.png');
% savefig('dischargeA.fig');
saveas(gcf,'C:/Users/maxim/OneDrive - The University of Sydney (Students)/THESIS/Thesis B figs/capbankdischarge.png');


% Estimated thrust from current level
Imax = max(current);
disp(['Peak current achieved (kA) = ' num2str(Imax / 1000)]);
Fmax = ((mu * Imax^2) / (4 * pi)) * (log(ra / rc) + 0.75);
disp(['Maximum theoretical thrust (N) = ' num2str(Fmax)]);

Ftheo = zeros(size(current));
for i = 1:length(current)
    Ftheo(i) = ((mu * current(i)^2) / (4 * pi)) * (log(ra / rc) + 0.75);
end

% Plot the force vs time
figure(2)
plot(t, Ftheo, 'Red', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Force (N)');
% set(gca, 'TickLabelInterpreter', 'latex');
% set(gca, 'TickLabelFormat', '%.1e');
grid on;
% saveas(gcf, 'pulseforce.png');
% savefig('pulseforce.fig');
saveas(gcf,'C:/Users/maxim/OneDrive - The University of Sydney (Students)/THESIS/Thesis B figs/Forcepulse.png');


% Function for ODE solver
function dudt = f(t, u)
    V = 1000;
    Ccap = 5 * (1 / (5 / (10e-6)));
    Rres = 1 / (5 / (5 * 0.1));
    Lind = 1 / (5 / (5 * (10e-9)));

    dudt = zeros(2, 1);
    dudt(1) = u(2);
    dudt(2) = (1 / Lind) * (V - (1 / Ccap) * u(1) - Rres * u(2));
end
