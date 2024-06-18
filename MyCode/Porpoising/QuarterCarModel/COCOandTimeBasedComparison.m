%% COCO and time based sim comparison
% In order for this script to work you first need to run the coco_cont
% script
format long
% Start by setting up the inputs
t = [0:0.01:10]';
% Our Inputs struct is already set up but want to use a value for a
% parameter that, according to COCO, will result in oscillations
paramtest = 'Kt';
% Define an array of test values
paramtestvalues = linspace(150000, 210000, 61);
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
    Inputs(Index4HB, 1) = paramtestvalues(i);
    % Time based test
    maxtime = 10;
    % Run the simulation and plot the results
    sol = ode45(@(t, x)SuspensionWithTime(t, x, Inputs), [0, maxtime], [-0.031; -0.008; 0; 0; 0; 0], odeset('RelTol', 1e-8)); % Simulate 
    DWFFloorResults = DWFFloor((Inputs(6, 1) + sol.y(1, :) + sol.y(2, :))', Inputs(9, 1), Inputs(10, 1), Inputs(11, 1))';
    % tiledlayout(1, 3); nexttile; plot(sol.x, sol.y(1, :)); xlabel('Time (secs)'); ylabel('Zs (m)');nexttile; plot(sol.x, sol.y(2, :)); ...
        % xlabel('Time (secs)');ylabel('Zu (m)');nexttile;plot(sol.x,Inputs(6, 1)+ sol.y(1, :) + sol.y(2, :)); xlabel('Time (secs)'); ylabel('Ride height (m)')
    
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
    [parampounique, ia, ib] = unique(vCarPO);
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
    TimePeriods = Po{1, 2}(2:end, find(strcmp(Po{1, 1}(1, :), 'po.period')));
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
plot(paramtestvalues, COCOzsAmplitudes)
xlabel(paramtest)
ylabel('Zs amplitude (m)')
legend('From time based sim', 'From COCO')
nexttile
plot(paramtestvalues, zuAmplitudes)
hold on
plot(paramtestvalues, COCOzuAmplitudes)
xlabel(paramtest)
ylabel('Zu amplitude (m)')
legend('From time based sim', 'From COCO')
nexttile
plot(paramtestvalues, rideheightAmplitudes)
hold on
plot(paramtestvalues, COCOrideheightAmplitudes)
xlabel(paramtest)
ylabel('Ride height amplitude (m)')
legend('From time based sim', 'From COCO')

figure 
plot(paramtestvalues, zsFreqs)
hold on
plot(paramtestvalues, zuFreqs)
plot(paramtestvalues, rideheightFreqs)
plot(paramtestvalues, COCOFreqs)
xlabel(paramtest)
ylabel('Frequency (Hz)')
legend('Zs from the time based sim', 'Zu from the time based sim', 'Ride height from the time based sim', 'From COCO')