%% Testing the qullllarter car model by running time series data through it
% Start by defining the time array and corresponding velocity array
% Want to accelerate from 0 to 320kph at 5ms^-2 to start
close all
clear all

vCar = [0:0.5:320]';
t = [0:0.1:(length(vCar)/10)-0.1]';

% Setting up the constant input terms for the suspension function
Inputs = [];
% Unsprung front mass (kg)
Inputs(1, 1) = 50;
% Unsprung rear mass (kg)
Inputs(2, 1) = 50;
% Sprung mass (kg)
Inputs(3, 1) = 659;
% Sprung mass moment of inertia (kgm^2)
Inputs(4, 1) = 710;
% Distance from front axle to CoG (m)
Inputs(5, 1) = 1.98;
% Distance from rear axle to CoG (m)
Inputs(6, 1) = 1.62;
% Front tyre vertical stiffness (N/m)
Inputs(7, 1) = 1.5 * 2.7 * 10^5;
% Front Suspension stiffness (N/m)
Inputs(8, 1) = 1.5 * 2 * 10^5;
% Rear tyre vertical stiffness (N/m)
Inputs(9, 1) =  1.5 * 0.9 * 2.7 * 10^5;
% Rear Suspension stiffness (N/m)
Inputs(10, 1) = 1.5 * 1.81 * 10^5;
% Front damping (Ns/m)
Inputs(11, 1) = 10434;
% Rear damping (Ns/m)
Inputs(12, 1) = 10992;
% Static front ride height (m)
Inputs(13, 1) = 0.0943;
% Static rear ride height (m)
Inputs(14, 1) = 0.11;
% Car speed (kph)
Inputs(15, 1) = 330;
% Front wing aero multiplier
Inputs(16, 1) = 0.254;
% Rear wing aero multiplier
Inputs(17, 1) = 0.636;

% Reading in the mat file that contains the parameters needed for the floor
% aero model
load("AeroValues.mat")

% Downforce from floor polynomial terms
Inputs(18, 1) = fittedDWFfunc.p00;
Inputs(19, 1) = fittedDWFfunc.p10;
Inputs(20, 1) = fittedDWFfunc.p01;
Inputs(21, 1) = fittedDWFfunc.p20;
Inputs(22, 1) = fittedDWFfunc.p11;
Inputs(23, 1) = fittedDWFfunc.p02;
Inputs(24, 1) = fittedDWFfunc.p30;
Inputs(25, 1) = fittedDWFfunc.p21;
Inputs(26, 1) = fittedDWFfunc.p12;
Inputs(27, 1) = fittedDWFfunc.p03;
Inputs(28, 1) = fittedDWFfunc.p31;
Inputs(29, 1) = fittedDWFfunc.p22;
Inputs(30, 1) = fittedDWFfunc.p13;
Inputs(31, 1) = fittedDWFfunc.p04;
Inputs(32, 1) = fittedDWFfunc.p32;
Inputs(33, 1) = fittedDWFfunc.p23;
Inputs(34, 1) = fittedDWFfunc.p14;
Inputs(35, 1) = fittedDWFfunc.p05;
Inputs(36, 1) = avgRHmean;
Inputs(37, 1) = avgRHstd;
Inputs(38, 1) = vCarmean;
Inputs(39, 1) = vCarstd;

sol = ode45(@(t, x)SuspensionWithTime(t, x, Inputs), [0, 5], [-0.015; -0.03; -0.06; -0.0054; 0; 0; 0; 0], odeset('RelTol', 1e-4)); % Simulate 

tiledlayout(3, 2); 
nexttile;
plot(sol.x, sol.y(1, :));
title('ZuF vs time');
fontsize(gca, 20, 'points')
nexttile;
plot(sol.x, sol.y(2, :));
title('ZuR vs time');
fontsize(gca, 20, 'points')
nexttile;
plot(sol.x, sol.y(3, :));
title('Zs vs time');
fontsize(gca, 20, 'points')
nexttile;
plot(sol.x, sol.y(4, :)); 
title('Theta vs time');
fontsize(gca, 20, 'points')
nexttile;
plot(sol.x, Inputs(13, 1) + sol.y(1, :) + sol.y(3, :) - Inputs(5, 1) .* tan(sol.y(4, :)));
title('hf vs time');
fontsize(gca, 20, 'points')
nexttile;
plot(sol.x, Inputs(14, 1) + sol.y(2, :) + sol.y(3, :) - Inputs(6, 1) .* tan(sol.y(4, :)));
title('hr vs time');
fontsize(gca, 20, 'points')


%% Testing stuff
% Now that the time based sim has been run can calculate the downforce
% generated by the floor, both at the front and the rear, for each point.
% This will be useful during the process of trying to create an aero
% function that works
% Start by defining the static floor geometry
% Defining the floor geometry
FloorLength = 3;
FloorIndices = 30001;

% Create an array of floor distances based on the floor length
FloorArray = linspace(0, FloorLength, FloorIndices)';

% Defining the x and y for the low point of the floor
Xlow = 1.7;
Ylow = Inputs(14, 1) - 0.03;

% Calculating the static rake angle
StaticRake = atan((Inputs(14, 1) - Inputs(13, 1)) / (Inputs(5, 1) + Inputs(6, 1)));

% Will form a cubic defining the floor profile
temp1 = zeros([4, 4]);
temp2 = zeros([4, 1]);

temp1(1, 4) = 1;
temp1(2, 1) = FloorLength^3;
temp1(2, 2) = FloorLength^2;
temp1(2, 3) = FloorLength;
temp1(2, 4) = 1;
temp1(3, 1) = Xlow^3;
temp1(3, 2) = Xlow^2;
temp1(3, 3) = Xlow;
temp1(3, 4) = 1;
temp1(4, 1) = 3 * Xlow^2;
temp1(4, 2) = 2 * Xlow;
temp1(4, 3) = 1;

temp2(1) = Inputs(13, 1) + 0.6 * tan(StaticRake);
temp2(2) = Inputs(14, 1);
temp2(3) = Ylow;
temp2(4) = 0;

% Solving
values = linsolve(temp1, temp2);
a = values(1);
b = values(2);
c = values(3);
d = values(4);

% Defining the floor profile
FloorProfile = a * FloorArray.^3 + b * FloorArray.^2 + c * FloorArray + d;

% Creating emtpy arrays to populate with values
DWFFloorF = zeros(length(sol.x), 1);
DWFFloorR = zeros(length(sol.x), 1);
TotalDWFFloor = zeros(length(sol.x), 1);
avgRH = zeros(length(sol.x), 1);
minfloorpoint = zeros(length(sol.x), 1);

% Running the for loop
for i = 1:length(sol.x)
    % Converting this static floor profile into the dynamic floor profile
    % First by applying the sprung mass displacment
    TempFloorProfile = FloorProfile + sol.y(3, i);
    
    % Then by applying the angle of the sprung mass
    THETAChange = sol.y(4, i) - atan((Inputs(14, 1) - Inputs(13, 1)) / (Inputs(5, 1) + Inputs(6, 1)));
    % Need to find the distance of each point in the floor from the point of
    % rotation to find it's veritcal displacment
    FloorDist2Rotation = FloorArray - (Inputs(5, 1) - 0.6);
    % Now finding the veritcal displacement due to the dynamic rotation of the
    % floor
    FloorVertChangeDue2Rotation = FloorDist2Rotation * tan(THETAChange);
    % Applying these vertical changes
    TempFloorProfile = TempFloorProfile + FloorVertChangeDue2Rotation;
    
    % Now that the dynamic floor profile has been found can compute the average
    % distance from the road to the floor
    avgRH(i) = mean(TempFloorProfile);
    minfloorpoint(i) = min(TempFloorProfile);
    
    % Normalising the vCar and average distance to road values 
    avgRHnormalised = (avgRH(i) - Inputs(36, 1)) ./ Inputs(37, 1);
    vCarnormalised = (Inputs(15, 1) - Inputs(38, 1)) ./ Inputs(39, 1);
    
    % Setting variables as x and y 
    z = avgRHnormalised;
    y = vCarnormalised;
    
    % Predicting the total downforce produced by the floor
    TotalDWFFloor(i) = Inputs(18, 1) + Inputs(19, 1).*z + Inputs(20, 1).*y + Inputs(21, 1).*z.^2 + Inputs(22, 1).*z.*y + ...
        + Inputs(23, 1).*y.^2 + Inputs(24, 1).*z.^3 + Inputs(25, 1).*z.^2.*y + Inputs(26, 1).*z.*y.^2 + Inputs(27, 1).*y.^3 ... 
        + Inputs(28, 1).*z.^3.*y + Inputs(29, 1).*z.^2.*y.^2 + Inputs(30, 1).*z.*y.^3 + Inputs(31, 1).*y.^4 ...
        + Inputs(32, 1).*z.^3.*y.^2 + Inputs(33, 1).*z.^2.*y.^3 + Inputs(34, 1).*z.*y.^4 + Inputs(35, 1).*y.^5;
    
    % Calculating the DWF on the front axle from the floor using the linear
    % approximation
    hF = Inputs(13, 1) + sol.y(3, i) + sol.y(1, i) - Inputs(5, 1) .* sol.y(4, i);
    DWFFloorF(i) = 0.45 * TotalDWFFloor(i);
    
    % From this obtaining the DWF on the rear axle from the floor
    DWFFloorR(i) = TotalDWFFloor(i) - DWFFloorF(i);
end

figure
tiledlayout(3, 2)
nexttile
plot(sol.x, DWFFloorF)
title('DWFFloorF')
nexttile
plot(sol.x, DWFFloorR)
title('DWFFloorR')
nexttile([1, 2])
plot(sol.x, TotalDWFFloor)
title('Total floor downforce')
nexttile
plot(sol.x, avgRH)
title('Average ride height value')
nexttile
plot(sol.x, minfloorpoint)
title('Minimum floor height')