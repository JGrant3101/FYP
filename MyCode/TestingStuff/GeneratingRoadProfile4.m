close all
clear all
clc
% Defining the PSD
B = 4.8 * 10^(-7); % Constant in PSD

vCar = 300; % Car speed in kph
vCarms = (vCar * 10^3) / (60 * 60); % Car speed in m/s

% Want to define our road length
RoadLength = 1500;
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

RoadProfileBandpass = bandpass(RoadProfile, [2, 15], fs);

figure
plot(TimePoints, RoadProfileBandpass)