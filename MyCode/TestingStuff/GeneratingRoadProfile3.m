% Define road length
close all
clear all
clc
RoadLength = 1500.01;
% Create evenly spaced points for road array
RoadPoints = linspace(0, RoadLength, RoadLength*100+1)';

% Define Car speed
vCar = 300;

% Convert road array to points in time
TimePoints = RoadPoints / ((vCar * 10^3)/(60*60));

% Now have a road profile defined as a function of time
% Find the sampling frequency of the road points
SamplingFreq = 1 / (TimePoints(2) - TimePoints(1));

N = length(TimePoints);

% Confirm total vertical distnace per metre of road
for i = 1:100
    RoadProfile = normrnd(0, 0.0133, N, 1);
    
    % plot(RoadPoints, RoadProfile)
    
    Differences = diff(RoadProfile);
    
    TotalVerticalMove(i, 1) = sum(abs(Differences));
    VerticalMovepermetre(i, 1) = TotalVerticalMove(i, 1) / RoadLength;
end


disp(mean(VerticalMovepermetre))
figure
plot(RoadPoints, RoadProfile)

