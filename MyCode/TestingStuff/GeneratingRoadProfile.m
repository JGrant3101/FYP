close all
clear all
clc

%% Attempt 1
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
    RoadProfile = normrnd(0, 50, N, 1);
    RoadProfile = bandpass(RoadProfile, [2,15], SamplingFreq);
    
    % plot(RoadPoints, RoadProfile)
    
    Differences = diff(RoadProfile);
    
    TotalVerticalMove(i, 1) = sum(abs(Differences));
    VerticalMovepermetre(i, 1) = TotalVerticalMove(i, 1) / RoadLength;
end


disp(mean(VerticalMovepermetre))
figure
plot(RoadPoints, RoadProfile)

xdft = fft(RoadProfile);
xdft = xdft(1:N/2+1);
psdx = (1/(SamplingFreq * N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:SamplingFreq/N:SamplingFreq/2;

figure
plot(freq, pow2db(psdx))

%% Attempt 2
close all
clear all
clc
% Defining the PSD
B = 4.8 * 10^(-7);
k = 2.1;
s = linspace(0, 1, 10001)';

Sg = B * s.^(-k);

loglog(s, Sg)
    
%% Attempt 3
close all 
clear all
clc

% Defining the PSD
B = 0.1;
k = 2;
f = linspace(0, 10, 101)';
vCar = (300 * 10^3) / (60 * 60);

Sh = B * vCar^(k-1) * f.^(-k);

plot(f, Sh)