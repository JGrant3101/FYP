%% Testing the qullllarter car model by running time series data through it
% Start by defining the time array and corresponding velocity array
% Want to accelerate from 0 to 320kph at 5ms^-2 to start
close all
clear all
vCar = [0:0.5:320]';
t = [0:0.1:(length(vCar)/10)-0.1]';
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
Inputs(7, 1) = 310;
% Upper downforce elements multiplier
Inputs(8, 1) = 0.365;
% Mean for Inverse Gaussian distribution
Inputs(9, 1) = 0.0001;
% Shape factor for Inverse Gaussian distribution
Inputs(10, 1) = 2.4;
% Scaling applied to Inverse Gaussian distribution
Inputs(11, 1) = 0.31*(500/9)^2;

maxtime = 10;

sol = ode45(@(t, x)SuspensionWithTime(t, x, Inputs), [0, maxtime], [-0.07; -0.02; 0; 0; 0; 0], odeset('RelTol', 1e-8)); % Simulate 
DWFFloorResults = DWFFloor((Inputs(6, 1) + sol.y(1, :) + sol.y(2, :))', Inputs(9, 1), Inputs(10, 1), Inputs(11, 1))';
    tiledlayout(1, 3); nexttile; plot(sol.x, sol.y(1, :)); title('Zs vs time'); xlabel('Time (seconds)'); ylabel('Zs (m)');
    nexttile; plot(sol.x, sol.y(2, :)); title('Zu vs time'); xlabel('Time (seconds)'); ylabel('Zu (m)'); nexttile;
    plot(sol.x, Inputs(6, 1) + sol.y(1, :) + sol.y(2, :)); title('Ride height vs time'); xlabel('Time (seconds)'); ylabel('h (m)');

figure
fontsize(gca, 20, 'points')
hold on
plot(sol.x, Inputs(6, 1) + sol.y(1, :) + sol.y(2, :));
title(['Ride height vs time for vCar = ', num2str(Inputs(7, 1)), ' kph']); 
xlabel('Time (seconds)'); 
ylabel('h (m)');

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
zsAmplitudes = zsMaxAverage - zsMinAverage;
zuAmplitudes = zuMaxAverage - zuMinAverage;
rideheightAmplitudes = rideheightMaxAverage - rideheightMinAverage;

% Printing these values
% temp = sprintf('From time based sims the amplitude of the Zs signal is %d', zsAmplitude);
% temp = temp + " (m)";
% disp(temp)
% temp = sprintf('From time based sims the amplitude of the Zu signal is %d', zuAmplitude);
% temp = temp + " (m)";
% disp(temp)
% temp = sprintf('From time based sims the amplitude of the ride height signal is %d', rideheightAmplitude);
% temp = temp + " (m)";
% disp(temp)

% Now to calculate the frequencies, the number of full periods in the time
% period is the same as the length of the smaller array of mins or maxes
% for each signal
Numzs = min([length(zsMaxIndices) length(zsMinIndices)]);
Numzu = min([length(zuMaxIndices) length(zuMinIndices)]);
Numrideheight = min([length(rideheightMaxIndices) length(rideheightMinIndices)]);

% Frequncy is then the number of periods used for these calculations
% divided by the total time over which those periods occur
if length(zsMaxIndices) > length(zsMinIndices)
    zsFreqs = Numzs / (sol.x(zsMaxIndices(end)) - sol.x(zsMaxIndices(1)));
else
    zsFreqs = Numzs / (sol.x(zsMinIndices(end)) - sol.x(zsMinIndices(1)));
end

if length(zuMaxIndices) > length(zuMinIndices)
    zuFreqs = Numzu / (sol.x(zuMaxIndices(end)) - sol.x(zuMaxIndices(1)));
else
    zuFreqs = Numzu / (sol.x(zuMinIndices(end)) - sol.x(zuMinIndices(1)));
end

if length(rideheightMaxIndices) > length(rideheightMinIndices)
    rideheightFreqs = Numrideheight / (sol.x(rideheightMaxIndices(end)) - sol.x(rideheightMaxIndices(1)));
else
    rideheightFreqs = Numrideheight / (sol.x(rideheightMinIndices(end)) - sol.x(rideheightMinIndices(1)));
end
return
% Now that all inputs have been set up can run the function in a for loop
% (which is required as ride heights calculated at the previous time step will have an impact on those at the next)
% to calculate the ride heights are each speed value.
% First want to initialise an x vector of all zeros
% x = zeros(4, 1);
% % Then want to initialise arrays which will be populated with results
% Zs = zeros(length(t), 1);
% Zu = zeros(length(t), 1);
% dotZs = zeros(length(t), 1);
% dotZu = zeros(length(t), 1);
% RideHeight = zeros(length(t), 1);
% RideHeight(:, 1) = 120;
% DWFUpper = zeros(length(t), 1);
% DWFFloorValues = zeros(length(t), 1);
% for i = 2:length(t)
%     % Assinging vCar
%     Inputs(7, 1) = vCar(i, 1);
%     % Calling the suspension function
%     [rhs, VarOfInterest] = Suspension(x, Inputs);
%     % Saving variables of interest
%     Zs(i, 1) = Zs(i-1, 1) + 0.1 * rhs(1, 1); 
%     Zu(i, 1) = Zu(i-1, 1) + 0.1 * rhs(2, 1);
%     %dotZs(i, 1) = dotZs(i-1, 1) + 0.1 * rhs(3, 1); 
%     %dotZu(i, 1) = dotZu(i-1, 1) + 0.1 * rhs(4, 1);
%     DWFUpper(i, 1) = VarOfInterest(1, 1);
%     DWFFloorValues(i, 1) = VarOfInterest(2, 1);
%     RideHeight(i, 1) = VarOfInterest(3, 1);
% 
%     % Creating new x array
%     x(1, 1) = Zs(i, 1);
%     x(2, 1) = Zu(i, 1);
%     x(3, 1) = rhs(3, 1);
%     x(4, 1) = rhs(4, 1);
% end