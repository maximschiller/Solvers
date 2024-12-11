clc
clear
clear all
close all
clf reset


% reading each sheet seperately into matlab
[~,sheet_name]=xlsfinfo('MPDtests_filter.xlsx')
% for k=1:numel(sheet_name)
%   data{k}=xlsread('MPDtests_filter.xlsx',sheet_name{k})
% end
 
k = 14;  % sheet value number

% sheet name 4,6,7,8,9,10,12,14
data = xlsread('MPDtests_filter.xlsx', sheet_name{k});
rows = size(data,1);
cols = size(data,2);
ntests = (cols + 1)/7; % number of measurements for each test
% time = [];
% capV = [];
% angle = [];
% accelx = [];
% accely = [];
% accelz = [];
testcounter = 0;

% Getting data from values for each dataset
for i = 1:ntests
%     time = [time data(:,testcounter + 1)];
%     capV = [capV data(:,testcounter + 2)];
%     angle = [angle data(:,testcounter + 3)];
%     accelx = [accelx data(:,testcounter + 4)];
%     accely = [accely data(:,testcounter + 5)];
%     accelz = [accelz data(:,testcounter + 6)];

    % storing data into cells
    time{i} = data(:,testcounter + 1);
    capV{i} = data(:,testcounter + 2);
    angle{i} = data(:,testcounter + 3);
    accelx{i} = data(:,testcounter + 4);
    accely{i} = data(:,testcounter + 5);
    accelz{i} = data(:,testcounter + 6);

    % calculating the acceleration magnitude
%     accelmag{i} = sqrt(accelx{i}.^2 + accely{i}.^2 + accelz{i}.^2);
    accelmag{i} = sqrt(accelx{i}.^2 + accely{i}.^2);

    testcounter = testcounter + 7;
end


for i = 1:ntests

    % removing nans
    for j = 1:length(accelmag{i})
        if isnan(accelmag{i}(j)) == 0 
                accelmagnan{i}(j) = accelmag{i}(j);
                timenan{i}(j) = time{i}(j);
                capVnan{i}(j) = capV{i}(j);
                anglenan{i}(j) = angle{i}(j);
                accelxnan{i}(j) = accelx{i}(j);
                accelynan{i}(j) = accely{i}(j);
                accelznan{i}(j) = accelz{i}(j);
        end
    end

    % removing zeros
    accelmagnan{i} = nonzeros(accelmagnan{i});
    timenan{i} = nonzeros(timenan{i});
    capVnan{i} = nonzeros(capVnan{i});
    anglenan{i} = nonzeros(anglenan{i});
    accelxnan{i} = nonzeros(accelxnan{i});
    accelynan{i} = nonzeros(accelynan{i});
    accelznan{i} = nonzeros(accelznan{i});


    % filtering the data
    accelmagnanfiltered{i} = smooth(accelmagnan{i},0.15,'loess');
    capVnanfiltered{i} = smooth(capVnan{i},0.15,'loess');
    anglenanfiltered{i} = smooth(anglenan{i},0.15,'loess');
    accelxnanfiltered{i} = smooth(accelxnan{i},0.15,'loess');
    accelynanfiltered{i} = smooth(accelynan{i},0.15,'loess');
    accelznanfiltered{i} = smooth(accelznan{i},0.15,'loess');
end


%% PLOTTING %%

figure(1)
for i = 1:ntests
    % Plotting
    plot(1:length(accelmagnan{i}),detrend(anglenanfiltered{i}), 'LineWidth',1.5)
    hold on
    L{i} = strcat('Test = ',num2str(i));
    legend(L, "Location","best")
end
grid on
grid minor
box on
xlabel("Counts");
ylabel("Angle (degrees)")
saveas(gcf,['C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\Thesis B figs\','MeasurementAngle_',sheet_name{k},'.png'])

figure(2)
for i = 1:ntests
    % Plotting
    plot(1:length(accelmagnan{i}),detrend(accelmagnanfiltered{i}), 'LineWidth',1.5)
%     plot(1:length(accelmagnan{i}),accelmagnanfiltered{i})
    hold on
    L{i} = strcat('Test = ',num2str(i));
    legend(L, "Location","best")
end
grid on
grid minor
box on
xlabel("Counts");
ylabel("Accelerometer Magnitude (g)")
saveas(gcf,['C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\Thesis B figs\','MeasurementAccelerometer_',sheet_name{k},'.png'])

figure(3)
for i = 1:ntests
    % Plotting
    plot(1:length(accelmagnan{i}),detrend(anglenan{i}), 'LineWidth',1.5)
    hold on
    L{i} = strcat('Test = ',num2str(i));
    legend(L, "Location","best")
end
grid on
grid minor
box on
xlabel("Counts");
ylabel("Angle (degrees)")
saveas(gcf,['C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\Thesis B figs\','MeasurementAngleunfiltered_',sheet_name{k},'.png'])

figure(4)
for i = 1:ntests
    % Plotting
    plot(1:length(accelmagnan{i}),detrend(accelmagnan{i}), 'LineWidth',1.5)
%     plot(1:length(accelmagnan{i}),accelmagnan{i})
    hold on
    L{i} = strcat('Test = ',num2str(i));
    legend(L, "Location","best")
end
grid on
grid minor
box on
xlabel("Counts");
ylabel("Accelerometer Magnitude (g)")
saveas(gcf,['C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\Thesis B figs\','MeasurementAccelerometerunfiltered_',sheet_name{k},'.png'])


%% Calculate thrust curves

mholder = 0.078;
mnozzles = [0.035 0.079 0.08 0.079]; % masses of the nozzles -> cylindrical, condi, flared, diverging
r = 0.26; % distance of thruster holder from centre of rotation
g = 9.81 % acceleration

% calculate the force

for i = 1:ntests
%     forcecyl{i} = accelmagnanfiltered{i}*mnozzles(1)*r;
%     forcecondi{i} = accelmagnanfiltered{i}*mnozzles(2)*r;
%     forceflared{i} = accelmagnanfiltered{i}*mnozzles(3)*r;
%     forcediverg{i} = accelmagnanfiltered{i}*mnozzles(4)*r;

    forcecyl{i} = accelmagnanfiltered{i}*mnozzles(1)*g;
    forcecondi{i} = accelmagnanfiltered{i}*mnozzles(2)*g;
    forceflared{i} = accelmagnanfiltered{i}*mnozzles(3)*g;
    forcediverg{i} = accelmagnanfiltered{i}*mnozzles(4)*g;
end

if k==12 || k==14
figure(5)
for i = 1:ntests
    % Plotting
    plot(1:length(accelmagnan{i}),detrend(forcecyl{i}), 'LineWidth',1.5)
    hold on
    L{i} = strcat('Test = ',num2str(i));
    legend(L, "Location","best")
end
grid on
grid minor
box on
xlabel("Counts");
ylabel("Force (N)")
saveas(gcf,['C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\Thesis B figs\','ForceCylinder_',sheet_name{k},'.png'])
end

if k == 4 || k==6
figure(6)
for i = 1:ntests
    % Plotting
    plot(1:length(accelmagnan{i}),detrend(forcecondi{i}), 'LineWidth',1.5)
    hold on
    L{i} = strcat('Test = ',num2str(i));
    legend(L, "Location","best")
end
grid on
grid minor
box on
xlabel("Counts");
ylabel("Force (N)")
saveas(gcf,['C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\Thesis B figs\','ForceCondi_',sheet_name{k},'.png'])
end

if k == 9 || k==10
figure(7)
for i = 1:ntests
    % Plotting
    plot(1:length(accelmagnan{i}),detrend(forceflared{i}), 'LineWidth',1.5)
    hold on
    L{i} = strcat('Test = ',num2str(i));
    legend(L, "Location","best")
end
grid on
grid minor
box on
xlabel("Counts");
ylabel("Force (N)")
saveas(gcf,['C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\Thesis B figs\','ForceFlared_',sheet_name{k},'.png'])
end

if k==7 || k==8
figure(8)
for i = 1:ntests
    % Plotting
    plot(1:length(accelmagnan{i}),detrend(forcediverg{i}), 'LineWidth',1.5)
    hold on
    L{i} = strcat('Test = ',num2str(i));
    legend(L, "Location","best")
end
grid on
grid minor
box on
xlabel("Counts");
ylabel("Force (N)")
saveas(gcf,['C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\Thesis B figs\','ForceDiverging_',sheet_name{k},'.png'])
end


%%  Difference in thrust and T/W Calcs

mthruster = 0.078;
mpropellant = 0.061;
mattachment = 0.068;
mbar = 0.079;
mcounterweight = 0.219;
mthyristor = 0.141;

hthruster = 0.026;
hcounterweight = 0.02;

% compute using thrust balance for condi values
if k == 3 || k == 5
    for i = 1:ntests
        hallThrustcondi{i} = ((mthruster + mpropellant + mattachment + mnozzles(2)) + mbar - (mcounterweight+mthyristor)*(hcounterweight/hthruster))*9.81*sind(detrend(anglenanfiltered{i}));
    end

    figure(9)
    for i = 1:ntests
        % Plotting
        plot(1:length(accelmagnan{i}),detrend(hallThrustcondi{i}), 'LineWidth',1.5)
        hold on
        L{i} = strcat('Test = ',num2str(i));
        legend(L, "Location","best")
    end
    grid on
    grid minor
    box on
    xlabel("Counts");
    ylabel("Force (N)")
    saveas(gcf,['C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\Thesis B figs\','hallforcecondi_',sheet_name{k},'.png'])
  
    maxthrustanglecondi{i} = max(detrend(hallThrustcondi{i}));
end

% Energy storage in capacitor bank
Vcap = 400;
Ccap = 10* 10^-6;
numcaps = 25;
E = Ccap*1350^2;

% getting the max value from each dataset in order to determine the maximum
% thrust produced
for i = 1:ntests

    forcecyl{i} = detrend(forcecyl{i});
    forcecondi{i} = detrend(forcecondi{i});
    forceflared{i} = detrend(forceflared{i});
    forcediverg{i} = detrend(forcediverg{i});

    maxthrustcyl{i} = max(forcecyl{i});
    maxthrustcondi{i} = max(forcecondi{i});
    maxthrustflared{i} = max(forceflared{i});
    maxthrustdiverg{i} = max(forcediverg{i});


    % computing T/W

    cylT_W{i} = maxthrustcyl{i}/mnozzles(1);
    condiT_W{i} = maxthrustcondi{i}/mnozzles(2);
    flaredT_W{i} = maxthrustflared{i}/mnozzles(3);
    divergT_W{i} = maxthrustdiverg{i}/mnozzles(4);

    % compute T/P
    cylT_P{i} = maxthrustcyl{i}/E;
    condiT_P{i} = maxthrustcondi{i}/E;
    flaredT_P{i} = maxthrustflared{i}/E;
    divergT_P{i} = maxthrustdiverg{i}/E;

end 

% concatenating data for outputting 
if k==12 || k==14
%     datacylthrust = mean(cell2mat(maxthrustcyl));
%     datacylT_W = mean(cell2mat(cylT_W));
%     datacylP_W = mean(cell2mat(cylT_P));
    datacylthrust = maxthrustcyl{1};
    datacylT_W = cylT_W{1};
    datacylP_W = cylT_P{1};
    data = [datacylthrust datacylT_W datacylP_W];
end

if k == 4 || k==6
%     datacondithrust = mean(cell2mat(maxthrustcondi));
%     datacondiT_W = mean(cell2mat(condiT_W));
%     datacondiP_W = mean(cell2mat(condiT_P));
    datacondithrust = maxthrustcondi{1};
    datacondiT_W = condiT_W{1};
    datacondiP_W = condiT_P{1};
    data = [datacondithrust datacondiT_W datacondiP_W];
end

if k == 9 || k==10
%     dataflaredthrust = mean(cell2mat(maxthrustflared));
%     dataflaredT_W = mean(cell2mat(flaredT_W));
%     dataflaredP_W = mean(cell2mat(flaredT_P));
    dataflaredthrust = maxthrustflared{1};
    dataflaredT_W = flaredT_W{1};
    dataflaredP_W = flaredT_P{1};
    data = [dataflaredthrust dataflaredT_W dataflaredP_W];
end

if k==7 || k==8
%     datadivergthrust = mean(cell2mat(maxthrustdiverg));
%     datadivergT_W = mean(cell2mat(divergT_W));
%     datadivergP_W = mean(cell2mat(divergT_P));
    datadivergthrust = maxthrustdiverg{1};
    datadivergT_W = divergT_W{1};
    datadivergP_W = divergT_P{1};
    data = [datadivergthrust datadivergT_W datadivergP_W];
end

d1 = 1.5/1000;
A1 = pi()*(d1/2)^2;
d2 = 4.7/1000;
A2 = pi()*(d2/2)^2;
rho = 1.293;


% calculating the mass flow rate through tube
% pressure difference 
dp = 101325 - 0.046;


v = sqrt(dp/(0.5*rho));

% % flow coefficient
% SGair = 1.03 * 10^-3;
% rho = 1.293;
% 
% Q = sqrt((2*dp)/rho)*(A2/sqrt(1 - (A2/A1)^2));
% Cf = Q*sqrt(SGair/dp);

Q = pi()*(d1/2)^2*v;

mdot = Q*rho;

% calculating specific impulse, impulse bit and efficiency
if k== 12 || k==14
    Isp = datacylthrust/(mdot*9.81);
    data = [data ; Isp 0 0];
end

if k == 4 || k==6
    Isp = datacondithrust/(mdot*9.81);
    data = [data ; Isp 0 0];
end

if k == 9 || k==10
    Isp = dataflaredthrust/(mdot*9.81);
    data = [data ; Isp 0 0];
end

if k == 7 || k==8
    Isp = datadivergthrust/(mdot*9.81);
    data = [data ; Isp 0 0];
end


%% calculating error in calibration of hall sensor

for i = 1:ntests
    angleerror{i} = detrend(anglenanfiltered{i}(1:100));
    meanangle(i) = mean(angleerror{i});
    stderror(i) = std(angleerror{i}) / sqrt(length(angleerror{i}));

end

data = [data ; mean(meanangle(i)) mean(stderror(i)) 0];