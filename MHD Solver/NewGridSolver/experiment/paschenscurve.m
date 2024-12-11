clc
clear

%% Paschens curve analysis
A = 112.5; % cm-1Torr -1
B = 2737.5; % Vcm-1 Torr-1
gamma =0.45;

% Breakdown voltage
pd = 0.567 % minimum breakdown voltage of air

% At p = 10^-4 pascals
p = 10;

d = [0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75]*10^-3;
% d = 7.5*10^-6;
for i = 1:length(d)
    Vb(i) = (B*p*d(i))/(log(A*p*d(i)) - log(log(1 + 1/gamma)));
end
plot(d,Vb)
% 
% 
% % Calculating values based on minimum breakdown voltage
% for i = 1:length(d)
%     p(i) = 0.567/d(i);
% end
% figure(2)
% plot(d,p)



%% Experimental breakdown of paschen law

d = linspace(0.001,1,100);
p = 0.01;
T = 273;

for i = 1:length(d)
    V(i) = 24.22*(293*p*d(i))/(760*T) + 6.08*((293*p*d(i))/(760*T))^0.5;
end
figure(3)
plot(d,V);