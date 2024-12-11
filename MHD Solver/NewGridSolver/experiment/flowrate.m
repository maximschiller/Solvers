clc
clear

d1 = 3/1000;
A1 = pi()*(d1/2)^2;
d2 = 2/1000;
A2 = pi()*(d2/2)^2;
rho = 1.293;


% pressure difference 
dp = 101325 - 0.046;


v = sqrt(dp/(0.5*rho));

% flow coefficient
SGair = 1.03 * 10^-3;
rho = 1.293;

Q = sqrt((2*dp)/rho)*(A2/sqrt(1 - (A2/A1)^2));
Cf = Q*sqrt(SGair/dp);

% Q = pi()*(d2/2)^2*v;
% 
% mdot = Q*rho;