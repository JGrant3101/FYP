% New script for doing the COCO vs time based so everything is just in one
% script 
clear all
close all
clc

%% COCO stuff

% Setting up the constant input terms for the suspension function
Inputs = [];
% Sprung mass (kg)
Inputs(1, 1) = 180;
% Unsprung mass (kg)
Inputs(2, 1) = 50;
% Suspension stiffness (N/m)
Inputs(3, 1) = 0.9*10^5;
% Suspension damping (Ns/m)
Inputs(4, 1) = 3775;
% Tyre vertical stiffness (N/m)
Inputs(5, 1) = 0.9*2.7 * 10^5;
% Static ride height (m)
Inputs(6, 1) = 0.1;
% vCar (kph)
Inputs(7, 1) = 280;
% Upper downforce elements multiplier
Inputs(8, 1) = 0.365;
% Mean for Inverse Gaussian distribution
Inputs(9, 1) = 0.0001;
% Shape factor for Inverse Gaussian distribution
Inputs(10, 1) = 2.4;
% Scaling applied to Inverse Gaussian distribution
Inputs(11, 1) = 0.31*(500/9)^2;

% Defining our initial x array
x0 = [-0.0371;-0.01;0;0];

% Defining the names of the inputs already written above
pnames = {'Ms' 'Mu' 'Ks' 'Cs' 'Kt' 'H' 'vCar' 'A' 'mew' 'lamda' 'scaling'};

% Initialising the COCO problem
prob = coco_prob();
prob = coco_set(prob, 'ep', 'NSA', true);

% Defining the functions
SystemSetup = {@Suspension}; %, @Suspension_dx, @Suspension_dp};

% Defining the function arguments
args = {SystemSetup{:}, x0, pnames, Inputs};

% Change the maximum step size
prob = coco_set(prob, 'cont', 'h_max', 1);

prob = coco_set(prob, 'cont', 'h_min', 1e-10);

% Increase the number of continuation steps
prob = coco_set(prob, 'cont', 'PtMX', 500);

% Calling function to vary a parameter and run a continuation
param = 'vCar';
Index = 7;

bdread0 = varyingparameters(param, [50 350], prob, args, Inputs);

% Plotting the HB point of the initial system
figure; clf; hold on
thm1 = struct('special', {{'HB', 'EP'}});
coco_plot_bd(thm1, 'Initial', 'vCar', '||x||_2')
title('Graph of the equilibrium point of the system as car speed varies')
grid on

% Further inspection of these HB points

labs = coco_bd_labs(bdread0, 'HB');

HBbd = cell(1, length(labs));
Po = cell(1, length(labs));

for i = 1:length(labs)
    prob2 = coco_prob();
    prob2 = coco_set(prob2, 'coll', 'NTST', 1140);
    prob2 = coco_set(prob2, 'po', 'bifus', true);
    prob2 = ode_HB2po(prob2, '', 'Initial', labs(i));
    prob2 = coco_set(prob2, 'cont', 'PtMX', 100);
    prob2 = coco_set(prob2, 'cont', 'ItMX', 2000);
    prob2 = coco_set(prob2, 'cont', 'h_min', 1e-10);
    prob2 = coco_set(prob2, 'cont', 'NAdapt', 1, 'h_max', 1);
    prob2 = coco_set(prob2, 'cont', 'TOL', 1e-3);
    prob2 = coco_set(prob2, 'coll', 'TOL', 1e-3);
    prob2 = coco_set(prob2, 'coco', 'TOL', 1e-3);
    prob2 = coco_set(prob2, 'EP', 'TOL', 1e-3);
    prob2 = coco_set(prob2, 'PO', 'TOL', 1e-3);
    prob2 = coco_set(prob2, 'ATLAS', 'TOL', 1e-3);

    % Tell COCO to store extra information about the periodic orbits
    prob2 = coco_add_slot(prob2, 'slot_bd_min_max', @slot_bd_min_max, [], 'bddat');

    Po{i} = coco(prob2, sprintf('po_run%d', i), [], 1, {param 'po.period'}, [50 350]);

    % Plot Equilibrium results for Zs
    paramep = coco_bd_col(bdread0, param); % Get the f parameter
    xep = coco_bd_col(bdread0, 'x'); % Get the x state
    lambdaep = coco_bd_col(bdread0, 'eigs'); % Get the x state
    stabep = any(real(lambdaep)>0) ;

    % Plot PO results for Zs
    bd2 = coco_bd_read(sprintf('po_run%d', i));
    parampo = coco_bd_col(bd2, param); % Get the current parameter
    x_max = coco_bd_col(bd2, 'x_max')';   % Get the state vector
    x_min = coco_bd_col(bd2, 'x_min')';   % Get the state vector
    stabpo = coco_bd_col(bd2, 'po.test.USTAB') == 0; % Get the stability

    % Plotting Zs results
    % Plotting the ep values
    figure
    xep = coco_bd_col(bdread0, 'x')';
    zsep = xep(:, 1);
    vCarep = coco_bd_col(bdread0, 'vCar')';

    plot(vCarep, zsep, 'LineWidth', 3)
    hold on

    % Plotting the PO values
    vCarPO = coco_bd_col(Po{i}, param);

    zs_max = x_max(1, :);
    zs_min = x_min(1, :);

    plot(vCarPO(stabpo), zs_max(stabpo),'.k', 'MarkerSize', 12) ;
    plot(vCarPO(~stabpo), zs_max(~stabpo),'.m', 'MarkerSize', 12) ;
    plot(vCarPO(stabpo), zs_min(stabpo),'.k', 'MarkerSize', 12) ;
    plot(vCarPO(~stabpo), zs_min(~stabpo),'.m', 'MarkerSize', 12) ;

    xlabel('vCar (kph)')
    ylabel('Zs (m)')
    title('EP and PO values of Zs across vCar range')
    xlim([290 360])
    fontsize(gca, 40, 'points')

    legend('Ep values', 'Stable PO results', 'Unstable PO results')
    % Plotting the Zu results
    % Plotting the ep values
    figure
    zuep = xep(:, 2);

    plot(vCarep, zuep, 'LineWidth', 3)
    hold on

    % Plotting the PO values
    zu_max = x_max(2, :);
    zu_min = x_min(2, :);

    plot(vCarPO(stabpo), zu_max(stabpo),'.k', 'MarkerSize', 12) ;
    plot(vCarPO(~stabpo), zu_max(~stabpo),'.m', 'MarkerSize', 12) ;
    plot(vCarPO(stabpo), zu_min(stabpo),'.k', 'MarkerSize', 12) ;
    plot(vCarPO(~stabpo), zu_min(~stabpo),'.m', 'MarkerSize', 12) ;

    xlabel('vCar (kph)')
    ylabel('Zu (m)')
    title('EP and PO values of Zu across vCar range')
    xlim([290 360])
    fontsize(gca, 40, 'points')

    legend('Ep values', 'Stable PO results', 'Unstable PO results')

    % Ride height results
    RideHeightEP = Inputs(6, 1) + zsep + zuep;
    RideHeightPO_Max = Inputs(6, 1) + zs_max + zu_max;
    RideHeightPO_Min = Inputs(6, 1) + zs_min + zu_min;

    figure
    plot(vCarep, RideHeightEP)
    hold on

    plot(vCarPO(stabpo), RideHeightPO_Max(stabpo),'.k') ;
    plot(vCarPO(~stabpo), RideHeightPO_Max(~stabpo),'.m') ;
    plot(vCarPO(stabpo), RideHeightPO_Min(stabpo),'.k') ;
    plot(vCarPO(~stabpo), RideHeightPO_Min(~stabpo),'.m') ;

    xlabel('vCar (kph)')
    ylabel('Ride height (m)')
    title('EP and PO values of ride height across vCar range')
    xlim([290 360])
    fontsize(gca, 40, 'points')

    legend('Ep values', 'Stable PO results', 'Unstable PO results')
end

%% Non COCO stuff

format long
% Start by setting up the inputs
t = [0:0.1:100]';
% Our Inputs struct is already set up but want to use a value for a
% parameter that, according to COCO, will result in oscillations
paramtest = param;
% Define an array of test values
paramtestvalues = linspace(307, 350, 44);
N = length(paramtestvalues);

% Defining arrays of 0s to add values to
zsAmplitudes = zeros(N, 1);
zuAmplitudes = zeros(N, 1);
rideheightAmplitudes = zeros(N, 1);

zsFreqs = zeros(N, 1);
zuFreqs = zeros(N, 1);
rideheightFreqs = zeros(N, 1);

COCOzsAmplitudes = zeros(N, 1);
COCOzuAmplitudes = zeros(N, 1);
COCOrideheightAmplitudes = zeros(N, 1);

COCOFreqs = zeros(N, 1);
for i = 1:N
    Inputs(Index, 1) = paramtestvalues(i);
    % Time based test
    maxtime = 10;
    % Run the simulation and plot the results
    sol = ode45(@(t, x)SuspensionWithTime(t, x, Inputs), [0, maxtime], [-0.031; -0.008; 0; 0; 0; 0], odeset('RelTol', 1e-8)); % Simulate 
    DWFFloorResults = DWFFloor((Inputs(6, 1) + sol.y(1, :) + sol.y(2, :))', Inputs(9, 1), Inputs(10, 1), Inputs(11, 1))';
    
    % In order to find the frequency and amplitude of the signal will look at
    % the last 5 seconds of the signal, at this point the system should have
    % settled into ocsillations and this range should give enough complete
    % periods of the signal to get reasonable amplitude and frequency values
    % Start by differentiating the signals with respect to time
    StartTimeIndex = find(sol.x == interp1(sol.x, sol.x, maxtime - 5, 'nearest'));
    zsdot = diff(sol.y(1, StartTimeIndex:end));
    zudot = diff(sol.y(2, StartTimeIndex:end));
    rideheightdot = diff(Inputs(6, 1) + sol.y(1, StartTimeIndex:end) + sol.y(2, StartTimeIndex:end));
    
    % Now that these differentials and double differentials have been obtained
    % want to find all the max and min points
    zsMaxIndices = find(zsdot(1:end-1)>0 & zsdot(2:end) < 0);
    zuMaxIndices = find(zudot(1:end-1)>0 & zudot(2:end) < 0);
    rideheightMaxIndices = find(rideheightdot(1:end-1)>0 & rideheightdot(2:end) < 0);
    
    zsMinIndices = find(zsdot(1:end-1)<0 & zsdot(2:end) > 0);
    zuMinIndices = find(zudot(1:end-1)<0 & zudot(2:end) > 0);
    rideheightMinIndices = find(rideheightdot(1:end-1)<0 & rideheightdot(2:end) > 0);
    
    % If a signal has a minimum first then to capture only full waves, which
    % will make frequency calculations easier, then we need the minimum array
    % to have one more term than the maximum array. THis is vice versa for the
    % opposite case.
    if zsMaxIndices(1) < zsMinIndices(1)
        if length(zsMaxIndices) == length(zsMinIndices)
            zsMinIndices(end) = [];
        end
    elseif zsMaxIndices(1) > zsMinIndices(1)
        if length(zsMaxIndices) == length(zsMinIndices)
            zsMaxIndices(end) = [];
        end
    end
    
    if zuMaxIndices(1) < zuMinIndices(1)
        if length(zuMaxIndices) == length(zuMinIndices)
            zuMinIndices(end) = [];
        end
    elseif zuMaxIndices(1) > zuMinIndices(1)
        if length(zuMaxIndices) == length(zuMinIndices)
            zuMaxIndices(end) = [];
        end
    end
    
    if rideheightMaxIndices(1) < rideheightMinIndices(1)
        if length(rideheightMaxIndices) == length(rideheightMinIndices)
            rideheightMinIndices(end) = [];
        end
    elseif rideheightMaxIndices(1) > rideheightMinIndices(1)
        if length(rideheightMaxIndices) == length(rideheightMinIndices)
            rideheightMaxIndices(end) = [];
        end
    end
    
    % Now all we are only capturing full oscillations and so calculate the
    % amplitude and frequency of these oscillations
    % First need to convert each of our index values so that they index the
    % original values arrays correctly
    zsMaxIndices = zsMaxIndices + StartTimeIndex;
    zuMaxIndices = zuMaxIndices + StartTimeIndex;
    rideheightMaxIndices = rideheightMaxIndices + StartTimeIndex;
    zsMinIndices = zsMinIndices + StartTimeIndex;
    zuMinIndices = zuMinIndices + StartTimeIndex;
    rideheightMinIndices = rideheightMinIndices + StartTimeIndex;
    
    % Now can find the average of the max and min values
    zsMaxAverage = mean(sol.y(1, zsMaxIndices));
    zuMaxAverage = mean(sol.y(2, zuMaxIndices));
    rideheightMaxAverage = mean(Inputs(6, 1) + sol.y(1, rideheightMaxIndices) + sol.y(2, rideheightMaxIndices));
    zsMinAverage = mean(sol.y(1, zsMinIndices));
    zuMinAverage = mean(sol.y(2, zuMinIndices));
    rideheightMinAverage = mean(Inputs(6, 1) + sol.y(1, rideheightMinIndices) + sol.y(2, rideheightMinIndices));
    
    % Now can calculate amplitudes
    zsAmplitudes(i) = zsMaxAverage - zsMinAverage;
    zuAmplitudes(i) = zuMaxAverage - zuMinAverage;
    rideheightAmplitudes(i) = rideheightMaxAverage - rideheightMinAverage;
    
    % Now to calculate the frequencies, the number of full periods in the time
    % period is the same as the length of the smaller array of mins or maxes
    % for each signal
    Numzs = min([length(zsMaxIndices) length(zsMinIndices)]);
    Numzu = min([length(zuMaxIndices) length(zuMinIndices)]);
    Numrideheight = min([length(rideheightMaxIndices) length(rideheightMinIndices)]);
    
    % Frequncy is then the number of periods used for these calculations
    % divided by the total time over which those periods occur
    if length(zsMaxIndices) > length(zsMinIndices)
        zsFreqs(i) = Numzs / (sol.x(zsMaxIndices(end)) - sol.x(zsMaxIndices(1)));
    else
        zsFreqs(i) = Numzs / (sol.x(zsMinIndices(end)) - sol.x(zsMinIndices(1)));
    end
    
    if length(zuMaxIndices) > length(zuMinIndices)
        zuFreqs(i) = Numzu / (sol.x(zuMaxIndices(end)) - sol.x(zuMaxIndices(1)));
    else
        zuFreqs(i) = Numzu / (sol.x(zuMinIndices(end)) - sol.x(zuMinIndices(1)));
    end
    
    if length(rideheightMaxIndices) > length(rideheightMinIndices)
        rideheightFreqs(i) = Numrideheight / (sol.x(rideheightMaxIndices(end)) - sol.x(rideheightMaxIndices(1)));
    else
        rideheightFreqs(i) = Numrideheight / (sol.x(rideheightMinIndices(end)) - sol.x(rideheightMinIndices(1)));
    end

    % Values from COCO
    % COCO has already calculated everything we need we just need to extract it
    % and report it to compare with the time based sim values
    % Want to start by interpolating the min and max values for Zs and Zu,
    % firstly need to ensure everything is only unique values
    [parampounique ia ib] = unique(parampo);
    x_max_unique = x_max(:, ia);
    x_min_unique = x_min(:, ia);
    
    zsmax = interp1(parampounique, x_max_unique(1, :), paramtestvalues(i));
    zsmin = interp1(parampounique, x_min_unique(1, :), paramtestvalues(i));
    
    zumax = interp1(parampounique, x_max_unique(2, :), paramtestvalues(i));
    zumin = interp1(parampounique, x_min_unique(2, :), paramtestvalues(i));
    
    % Can now calculate amplituds and print these values
    COCOzsAmplitudes(i) = zsmax - zsmin;
    COCOzuAmplitudes(i) = zumax - zumin;
    % Will simply sum these together to get the overall ride heigh amplitude
    COCOrideheightAmplitudes(i) = COCOzsAmplitudes(i) + COCOzuAmplitudes(i);

    % Now want to extract the period from our COCO data in order to get the
    % frequency of the oscillations at this point according to COCO, this will
    % just be one value
    TimePeriods = Po{1, 1}(2:end, find(strcmp(Po{1, 1}(1, :), 'po.period')));
    TimePeriods = cell2mat(TimePeriods);
    TimePeriods_unique = TimePeriods(ia);
    TimePeriod = interp1(parampounique, TimePeriods_unique, paramtestvalues(i));
    
    % Can now get the frequency and print it
    COCOFreqs(i) = 1/TimePeriod;
end

% Plotting the results
figure
tiledlayout(1, 3)
nexttile
plot(paramtestvalues, zsAmplitudes)
hold on
fontsize(gca, 20, 'points')
plot(paramtestvalues, COCOzsAmplitudes)
xlabel([paramtest, ' (kph)'])
ylabel('Zs amplitude (m)')
legend('From time based sim', 'From COCO')
nexttile
plot(paramtestvalues, zuAmplitudes)
hold on
plot(paramtestvalues, COCOzuAmplitudes)
xlabel([paramtest, ' (kph)'])
fontsize(gca, 20, 'points')
ylabel('Zu amplitude (m)')
legend('From time based sim', 'From COCO')
nexttile
plot(paramtestvalues, rideheightAmplitudes)
hold on
plot(paramtestvalues, COCOrideheightAmplitudes)
xlabel([paramtest, ' (kph)'])
ylabel('Ride height amplitude (m)')
legend('From time based sim', 'From COCO')
fontsize(gca, 20, 'points')

figure 
plot(paramtestvalues, zsFreqs)
hold on
plot(paramtestvalues, zuFreqs)
plot(paramtestvalues, rideheightFreqs)
plot(paramtestvalues, COCOFreqs)
xlabel([paramtest, ' (kph)'])
ylabel('Frequency (Hz)')
legend('Zs from the time based sim', 'Zu from the time based sim', 'Ride height from the time based sim', 'From COCO')
title('Comparison of oscillation frequencies for time based simulations and COCO')
fontsize(gca, 20, 'points')

figure
plot(paramtestvalues, rideheightAmplitudes)
hold on
plot(paramtestvalues, COCOrideheightAmplitudes)
xlabel([paramtest, ' (kph)'])
ylabel('Ride height amplitude (m)')
legend('From time based sim', 'From COCO')
title('Comparison of oscillation amplitudes for time based simulations and COCO')
fontsize(gca, 20, 'points')

%% Calculations based on the results
% For every value caclulate the absolute difference between the time based
% sim and COCO output value
zsAmplitudesDiff = COCOzsAmplitudes - zsAmplitudes;
zuAmplitudesDiff = COCOzuAmplitudes - zuAmplitudes;
rideheightAmplitudesDiff = COCOrideheightAmplitudes - rideheightAmplitudes;

rideheightFreqsDiff = COCOFreqs - rideheightFreqs;

% Recalculating these as percentage differences
zsAmplitudesPercentDiff = (zsAmplitudesDiff ./ zsAmplitudes) * 100;
zuAmplitudesPercentDiff = (zuAmplitudesDiff ./ zuAmplitudes) * 100;
rideheightAmplitudesPercentDiff = (rideheightAmplitudesDiff ./ rideheightAmplitudes) * 100;
rideheightFreqsPercentDiff = (rideheightFreqsDiff ./ rideheightFreqs) * 100;

% Finding the average of these percentage differences
avgzsAmplitudesPercentDiff = mean(zsAmplitudesPercentDiff);
avgzuAmplitudesPercentDiff = mean(zuAmplitudesPercentDiff);
avgrideheightAmplitudesPercentDiff = mean(rideheightAmplitudesPercentDiff);
avgrideheightFreqsPercentDiff = mean(rideheightFreqsPercentDiff);

%% Section for functions
function bdread0 = varyingparameters(param2vary, range, prob, args, Inputs)
    variableargs = {1, param2vary, range};

    bd0 = coco(prob, 'Initial', @ode_isol2ep, args{:}, variableargs{:});
    
    
    bdread0 = coco_bd_read('Initial'); % Get the Omega parameter
    [~,column] = find(strcmp(bdread0,param2vary));
    IndependentVar_val = cell2mat([bdread0(2:end,column)])' ;

    x_val = cell2mat([bdread0(2:end,28)]') ;
    figure ; plot(IndependentVar_val, x_val(1:2,:)) ;
    hep = Inputs(6, 1) + x_val(1,:) + x_val(2,:);
    hold on
    plot(IndependentVar_val, hep)
    xlabel(param2vary)
    ylabel('Equilibrium values (m)')
    legend('Zs (m)', 'Zu (ms)', 'Ride height (m)')
end
