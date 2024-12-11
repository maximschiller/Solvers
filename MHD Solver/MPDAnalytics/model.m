clc
clear


% Creating dataset
ra = 0.01:0.01:0.5;
rc = 0.01:0.01:0.5;
J = 0:500:10000;
[ra, rc, J] = meshgrid(ra, rc, J);

mu0 = 4*pi*10^-7;
% Thrust Maecker model
Tsf = (mu0/(4*pi))*(log(ra./rc) + 3/4).*J.^2;

% Creating plot
% ax = fig.add_subplot(1,1,1,'projection','3d');
scatter3(ra(:), rc(:), J(:), 20, Tsf(:), 'filled', 'MarkerFaceAlpha', .8);
colorbar();
% Set axes label
xlabel('ra', 'labelpad', 5);
ylabel('rc', 'labelpad', 5);
zlabel('T/Jsq', 'labelpad', 5);