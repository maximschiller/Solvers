clear all; close all; clc;
format long;

PI = 2*asin(1);
XMAX = 300;
YMAX = 100;

IMAX = 101; % No. of nodes in ksi direction i.e. (IMAX-1 Div)
JMAX = 26; % Initial No. of nodes in eta direction i.e. (JMAX-1 Div)

% Constants and arrays initialization
thickness = 0.12;

% Memory allocation and initialization
for i = 1:(IMAX+2)
    for j = 1:IMAX
        x(i,j) = 0.001;
        y(i,j) = 0.001;
        P(i,j) = 0.0;
        Q(i,j) = 0.0;
    end
end

% Main Program...Lower half Boundary Points
r = 0.8; % Ratio in ksi direction in the cut
h = 3*(1-r)/(1-power(r,(IMAX-1)/5)); % smallest h in ksi direction in the cut
dh = 0;
for i = 1:(IMAX-1)/2
    if i < (IMAX-1)/5
        x(i,1) = 4-dh;
        y(i,1) = 0;
        x(i,JMAX) = 4-dh;
        y(i,JMAX) = -2;
        dh = dh + h*power(r, i);
    elseif i < 7*(IMAX-1)/20
        x(i,1) = 0.5 + 0.5*cos(0.5*PI*(i-(IMAX-1)/5)/(7*(IMAX-1)/20-(IMAX-1)/5));
        y(i,1) = -5*thickness*(0.2948*power(x(i,1),0.5)-0.126*x(i,1)-0.3516*power(x(i,1),2)+0.2843*power(x(i,1),3)-0.1015*power(x(i,1),4));
        x(i,JMAX) = 1 - 3*sin((i-(IMAX-1)/5)*PI/(2*((IMAX-1)/2-(IMAX-1)/5)))^1.3;
        y(i,JMAX) = -2*cos((i-(IMAX-1)/5)*PI/(2*((IMAX-1)/2-(IMAX-1)/5)));
    else
        x(i,1) = 0.5-0.5*sin(0.5*PI*(i-(7*(IMAX-1)/20))/(7*(IMAX-1)/20-(IMAX-1)/5));
        y(i,1) = -5*thickness*(0.2948*power(x(i,1),0.5)-0.126*x(i,1)-0.3516*power(x(i,1),2)+0.2843*power(x(i,1),3)-0.1015*power(x(i,1),4));
        x(i,JMAX) = 1 - 3*sin((i-(IMAX-1)/5)*PI/(2*((IMAX-1)/2-(IMAX-1)/5)))^1.3;
        y(i,JMAX) = -2*cos((i-(IMAX-1)/5)*PI/(2*((IMAX-1)/2-(IMAX-1)/5)));
    end
end

%% test

clc
clear

imax = 101;
r =0.8;
h = (3*(1-r))/(1-r^((imax-1)/5));

dh = 0;
x = [];
k = [];
for i = 0:(imax-1)/2
    if i<(imax-1)/5
        x = [x 4 - dh];
        k = [k 1];
        dh = dh + h*r^i;
    elseif i < 7*(imax-1)/20
        x = [x 0.5 + 0.5*cos(0.5*pi()*(i-(imax-1)/5)/(7*(imax-1)/20-(imax-1)/5))];
    else 
        x = [x 0.5-0.5*sin(0.5*pi()*(i-(7*(imax-1)/20))/(7*(imax-1)/20-(imax-1)/5))];
    end
end
