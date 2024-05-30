% Running simulations with a not perfectly flat road input
%% Running the simulations
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
Inputs(4, 1) = 3400;
% Tyre vertical stiffness (N/m)
Inputs(5, 1) = 0.9*2.7 * 10^5;
% Static ride height (m)
Inputs(6, 1) = 0.1;
% vCar (kph)
Inputs(7, 1) = 250;
% Upper downforce elements multiplier
Inputs(8, 1) = 0.365;
% Mean for Inverse Gaussian distribution
Inputs(9, 1) = 0.0001;
% Shape factor for Inverse Gaussian distribution
Inputs(10, 1) = 2.4;
% Scaling applied to Inverse Gaussian distribution
Inputs(11, 1) = 0.31*(500/9)^2;

maxtime = 15;

% Defining the road profile
% Defining the PSD
B = 4.8 * 10^(-9); % Constant in PSD

vCarms = (Inputs(7, 1) * 10^3) / (60 * 60); % Car speed in m/s

% Defining the road input
% Want to define our road length
RoadLength = vCarms * (maxtime + 0.5);
% Splitting into road domain points
RoadPoints = 0:0.1:RoadLength;
RoadPoints = RoadPoints';
% Converting this into time points
TimePoints = RoadPoints / vCarms;
N = length(TimePoints);
fs = 1 / (TimePoints(2) - TimePoints(1));

% Defining an array of random normal values from distrib with E[0] Var[1]
randomnormals = normrnd(0, 1, N, 1);
DvbyDT = 2 * pi * sqrt(vCarms * B) * randomnormals;

RoadProfile = cumtrapz(TimePoints, DvbyDT);

figure
plot(TimePoints, RoadProfile)

% Filter out the very high frequency noise and low frequency macro changes
RoadProfileBandpass = bandpass(RoadProfile, [2, 15], fs);
% Centre around the first point
RoadProfileBandpass = RoadProfileBandpass - RoadProfileBandpass(1);

% Forming this into a 2d array lookup table
RoadProfileInput = [TimePoints, RoadProfileBandpass];

figure
plot(TimePoints, RoadProfileBandpass)

sol = ode45(@(t, x)SuspensionWithTime(t, x, Inputs, RoadProfileInput), [0, maxtime], [-0.031; -0.008; 0; 0; 0; 0], odeset('RelTol', 1e-8)); % Simulate

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