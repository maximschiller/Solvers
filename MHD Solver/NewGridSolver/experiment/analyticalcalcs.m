clc
clear

mu = (10^-7)*4*pi();

% Maecker Model
ra_rc = linspace(1,10,10);
current = linspace(0,20000,500);

T = [];
for i = 1:length(ra_rc)
    for j = 1:length(current)
        T(i,j) = (mu/(4*pi()))*(log(ra_rc(i)) + 0.75)*current(j)^2;
    end
    plot(current,T(i,:), 'LineWidth', 2);
    Legend{i}=strcat('ra/rc = ', num2str(ra_rc(i)));
    hold on
end
grid on
grid minor
box on 
legend(Legend, Location="best")
xlabel('Current (A)')
ylabel('Maecker Thrust (N)')
saveas(gcf,'C:/Users/maxim/OneDrive - The University of Sydney (Students)/THESIS/Thesis B figs/MaeckerThrust.png');


current = linspace(0,20000,10);
ra_r = linspace(1,10,500);
figure(7)
% variation in thrust with gap size
T = [];
for i = 1:length(current)
    for j = 1:length(ra_rc)
        T(i,j) = (mu/(4*pi()))*(log(ra_rc(j)) + 0.75)*current(i)^2;
    end
    plot(ra_rc,T(i,:), 'LineWidth', 2);
    Legend{i}=strcat('Current (A) = ', num2str(current(i)));
    hold on
end
grid on
grid minor
box on 
legend(Legend, Location="best")
xlabel('Gap Size Ratio (ra/rc)')
ylabel('Maecker Thrust (N)')
saveas(gcf,'C:/Users/maxim/OneDrive - The University of Sydney (Students)/THESIS/Thesis B figs/MaeckerThrust_gapsize.png');



%% CURRENT REQUIRED 

Frange = linspace(1*10^-3,100*10^-3,500);
ra = 10/1000
rc = 5/1000

for i = 1:length(Frange)
    I(i) = sqrt((Frange(i)*4*pi)/(mu*(log(ra/rc + 0.75))));
end

figure(3)
plot(Frange, I, 'LineWidth', 2)
grid on
grid minor
box on 
ylabel('Current Required (A)')
xlabel('Thrust (N)')
saveas(gcf,'C:/Users/maxim/OneDrive - The University of Sydney (Students)/THESIS/Thesis B figs/currentrequired.png');


%% BREAKDOWN VOLTAGE
% Paschens Law analysis

% Constants for air
a = 112.5; % kPa/cm
b = 2737.5; % V/kPa.cm
% a = 15; % /cm.Torr
% b = 365; % V/cm.Torr
gamma = 10^-2;
c = log(a / log(1 + 1/gamma));

% at various gap spacings
pressure = linspace(10^-5,10^6,1000)/1000;
% pressure = linspace(10^-1,10^3,1000); % Torr
ra = 20; % cm
d = [];
for i = 1:length(ra_rc)
    rc(i) = ra/ra_rc(i);
    d(i) = ra - rc(i);
end

figure(2)
% Breakdown voltage
for i = 2:length(ra_rc)
    for j = 1:length(pressure)
%         Vbreakdown(i,j) = b*pressure(j)*d(i) / (c + log(pressure(j)*d(i)));
        Vbreakdown(i,j) = (b*pressure(j)*d(i))/(log(a*pressure(j)*d(i)) - log(log(1 + 1/gamma)));
    end
    plot(pressure,Vbreakdown(i,:))
%     semilogx(pressure,Vbreakdown(i,:));
    hold on
end
grid on
grid minor
box on 
xlabel('Pressure (kPa)')
ylabel('Breakdown Voltage (V)')
saveas(gcf,'C:/Users/maxim/OneDrive - The University of Sydney (Students)/THESIS/Thesis B figs/Vbreakdown.png');

figure(4)
d = linspace(0.1,1,10);
for i = 1:length(d)
for j  = 1:length(pressure)
    Vb(i,j) = (b*pressure(j)*d(i))/(log(a*pressure(j)*d(i)) - log(log(1 + 1/gamma)));
end

plot(pressure,Vb(i,:))
hold on
end


grid on
grid minor
box on 
xlabel('Pressure (kPa)')
ylabel('Breakdown Voltage (V)')


% Breakdown field strength



% Minimum voltage
