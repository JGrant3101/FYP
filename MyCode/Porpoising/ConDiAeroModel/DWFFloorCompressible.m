% Creating a floor downforce function based on geometry and airflow as
% opposed to just a numerical approximation
clear all
close all
clc

%% Defining the floor geometry
% Going to start by redefining everything so I can just tweak this section
% Start by defining the floor length and width
FloorLength = 3;
FloorWidth = 1.5;
FloorIndices = 30001;

% Create an array of floor distances based on the floor length
FloorArray = linspace(0, FloorLength, FloorIndices)';

% Start by defining static front and rear ride heights
StaticRRH = 0.11;
FloorAngle = 0.25;
StaticFRH = StaticRRH - FloorLength * sind(FloorAngle);
%StaticFRH = 0.0943;


% Defining the x and y for the low point of the floor
Xlow = 1.7;
Ylow = 0.08;
% Will form a quadratic defining the floor profile
temp1 = zeros([4, 4]);
temp2 = zeros([4, 1]);

% Lowering the floor here while working on code
StaticFRH = StaticFRH - 0.074;
StaticRRH = StaticRRH - 0.074;
Ylow = Ylow - 0.074;

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

temp2(1) = StaticFRH;
temp2(2) = StaticRRH;
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

% While working on this code will plot this
figure
plot(FloorArray, FloorProfile)
xlim([0, FloorLength])
ylim([0, StaticRRH])
%axis equal

%% Doing the airflow stuff 
% Defining the inlet air conditions
Ti = 293.15; % This is in Kelvin so is 20 degrees celcius
Pi = 10^5; % Atmospheric pressure set to 1 bar
R = 287;
gamma = 1.4;
density = Pi / (287 * Ti);
% Will assume density is constant throughout
SoSi = sqrt(gamma * R * Ti);

% Define the car speed which is also inlet air speed
vCar = 320; % This is in kph will convert to m/s
vCar = (vCar * 10^3)/3600;

Mi = vCar / SoSi;

% Defining empty arrays to be populated
airspeeds = zeros([FloorIndices, 1]);
airMachs = zeros([FloorIndices, 1]);
airPressures = zeros([FloorIndices, 1]);
airTemps = zeros([FloorIndices, 1]);
Sos = zeros([FloorIndices, 1]);

airspeeds(1) = vCar;
airMachs(1) = Mi;
airPressures(1) = Pi;
airTemps(1) = Ti;
Sos(1) = SoSi;

% Defining an array of areas
A = FloorProfile * FloorWidth;

% Defining the constants to be used throughout
constant1 = Ti * (1 + 0.5 * (gamma - 1) * Mi^2);
constant2 = Ti / (Pi^((gamma - 1) / gamma));

% Pre setting the shock flag
shock = false;

% Setting the shock mach multiplier
ShockMulti = 2.1;

% Running the for loop across the length of the floor
for i = 2:FloorIndices
    if ~shock
        % Finding the new air speed
        deltaA = A(i) - A(i-1);
        deltaairspeeds = (airspeeds(i-1) * deltaA) / (A(i-1) * (airMachs(i-1)^2 - 1));
    
        airspeeds(i) = airspeeds(i-1) + deltaairspeeds;
        
        % Finding the new air temp
        airTemps(i) = constant1 - (((gamma - 1) * airspeeds(i)^2) / (2 * gamma * R));
    
        % Finding the new mach number
        Sos(i) = sqrt(gamma * R * airTemps(i));
        airMachs(i) = airspeeds(i) / Sos(i);
    
        % Once mach is found want to check if it is greater than 1 and if so
        % correct it as a normal shock would have occurred
        if airMachs(i) > 0.99
            shock = true;
            airMachs(i) = ShockMulti * airMachs(i);
        end

        % Finding the new air pressure
        airPressures(i) = (airTemps(i) / constant2) ^ (gamma / (gamma - 1));
    else
        % This part runs if mach number was above 1
        airMachs(i) = sqrt((1 + 0.5 * (gamma - 1) * airMachs(i-1)^2) / (gamma * airMachs(i-1)^2 - 0.5 * (gamma - 1)));

        airMachs(i) = airMachs(i) * 1/ShockMulti;
    
        airTemps(i) = (airTemps(i-1) * (1 + 0.5 * (gamma - 1) * airMachs(i-1)^2)) / (1 + 0.5 * (gamma - 1) * airMachs(i)^2);

        Sos(i) = sqrt(gamma * R * airTemps(i));

        airspeeds(i) = airMachs(i) * Sos(i);

        airPressures(i) = (airPressures(i-1) * (1 + gamma * airMachs(i-1)^2)) / (1 + gamma * airMachs(i)^2);

        shock = false;
    end
end

% Plotting the variation of area, air speed, mach number and pressure along
% the length of the floor
figure 
tiledlayout(2, 2)
nexttile
plot(FloorArray, A)
xlabel('Distance along the floor (m)')
ylabel('Floor cross sectional area (m^2)')

nexttile
plot(FloorArray, airspeeds)
xlabel('Distance along the floor (m)')
ylabel('Air speed under the floor (m/s)')

nexttile
plot(FloorArray, airMachs)
xlabel('Distance along the floor (m)')
ylabel('Mach number under the floor')

nexttile
plot(FloorArray, airPressures)
xlabel('Distance along the floor (m)')
ylabel('Air pressure under the floor (kg/m^2)')

% Finding the difference between the pressure under the floor and the
% ambient pressure
deltaP = Pi - airPressures;

% Integrating this over the floor area to find the overall downforce
% generated by the floor
Downforce = sum(deltaP * FloorWidth * (FloorArray(2) - FloorArray(1)));

% Converting this to kilograms of downforce
DownforceKG = Downforce / 9.81;

% Printing the output value
disp(DownforceKG)

%% Running the sweeps
% The above is giving downforce values that look vaguely reasonable if a
% bit high, want to now look into creating 3D plots of downforce values
% from various rear ride heights and car speeds, will not be considering
% the angle of the floor for now
% Want to define an array of car speeds and ride heights to test
vCararray = linspace(50, 350, 61)';
Diff2staticarray = linspace(0, -0.075, 76)';

% Creating an empty downforce array
DownforceVals = zeros(length(vCararray), length(Diff2staticarray));
avgdist2road = zeros(length(vCararray), length(Diff2staticarray));;

% Running a double for loop to populate the downforce array
for i = 1:length(vCararray)
    for j = 1:length(Diff2staticarray)
        % Start by defining static front and rear ride heights
        StaticRRH = 0.11;
        FloorAngle = 0.25;
        StaticFRH = StaticRRH - FloorLength * sind(FloorAngle);
        
        % Defining the x and y for the low point of the floor
        Xlow = 1.7;
        Ylow = 0.08;

        % Will form a quadratic defining the floor profile
        temp1 = zeros([4, 4]);
        temp2 = zeros([4, 1]);
        
        
        % Lowering the floor here while working on code
        StaticFRH = StaticFRH + Diff2staticarray(j);
        StaticRRH = StaticRRH + Diff2staticarray(j);
        Ylow = Ylow + Diff2staticarray(j);
        
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
        
        temp2(1) = StaticFRH;
        temp2(2) = StaticRRH;
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

        % Finding the average distance between the floor and the road
        avgdist2road(i, j) = mean(FloorProfile);

        % Defining the inlet air conditions
        Ti = 293.15; % This is in Kelvin so is 20 degrees celcius
        Pi = 10^5; % Atmospheric pressure set to 1 bar
        R = 287;
        gamma = 1.4;
        density = Pi / (287 * Ti);
        % Will assume density is constant throughout
        SoSi = sqrt(gamma * R * Ti);
        
        % Define the car speed which is also inlet air speed
        vCar = vCararray(i); % This is in kph will convert to m/s
        vCar = (vCar * 10^3)/3600;
        
        Mi = vCar / SoSi;
        
        % Defining empty arrays to be populated
        airspeeds = zeros([FloorIndices, 1]);
        airMachs = zeros([FloorIndices, 1]);
        airPressures = zeros([FloorIndices, 1]);
        airTemps = zeros([FloorIndices, 1]);
        Sos = zeros([FloorIndices, 1]);
        
        airspeeds(1) = vCar;
        airMachs(1) = Mi;
        airPressures(1) = Pi;
        airTemps(1) = Ti;
        Sos(1) = SoSi;
        
        % Defining an array of areas
        A = FloorProfile * FloorWidth;
        
        % Defining the constants to be used throughout
        constant1 = Ti * (1 + 0.5 * (gamma - 1) * Mi^2);
        constant2 = Ti / (Pi^((gamma - 1) / gamma));
        
        % Pre setting the shock flag
        shock = false;
        
        % Setting the shock mach multiplier
        ShockMulti = 2.5;
        
        % Running the for loop across the length of the floor
        for k = 2:FloorIndices
            if ~shock
                % Finding the new air speed
                deltaA = A(k) - A(k-1);
                deltaairspeeds = (airspeeds(k-1) * deltaA) / (A(k-1) * (airMachs(k-1)^2 - 1));
            
                airspeeds(k) = airspeeds(k-1) + deltaairspeeds;
                
                % Finding the new air temp
                airTemps(k) = constant1 - (((gamma - 1) * airspeeds(k)^2) / (2 * gamma * R));
            
                % Finding the new mach number
                Sos(k) = sqrt(gamma * R * airTemps(k));
                airMachs(k) = airspeeds(k) / Sos(k);
            
                % Once mach is found want to check if it is greater than 1 and if so
                % correct it as a normal shock would have occurred
                if airMachs(k) > 0.99
                    shock = true;
                    airMachs(k) = ShockMulti * airMachs(k);
                end
        
                % Finding the new air pressure
                airPressures(k) = (airTemps(k) / constant2) ^ (gamma / (gamma - 1));
            else
                % This part runs if mach number was above 1
                airMachs(k) = sqrt((1 + 0.5 * (gamma - 1) * airMachs(k-1)^2) / (gamma * airMachs(k-1)^2 - 0.5 * (gamma - 1)));
        
                airMachs(k) = airMachs(k) * 1/ShockMulti;
            
                airTemps(k) = (airTemps(k-1) * (1 + 0.5 * (gamma - 1) * airMachs(k-1)^2)) / (1 + 0.5 * (gamma - 1) * airMachs(k)^2);
        
                Sos(k) = sqrt(gamma * R * airTemps(k));
        
                airspeeds(k) = airMachs(k) * Sos(k);
        
                airPressures(k) = (airPressures(k-1) * (1 + gamma * airMachs(k-1)^2)) / (1 + gamma * airMachs(k)^2);
        
                shock = false;
            end
        end
        
        % Finding the difference between the pressure under the floor and the
        % ambient pressure
        deltaP = Pi - airPressures;

        % Integrating this over the floor area to find the overall downforce
        % generated by the floor
        DownforceVals(i, j) = sum(deltaP * FloorWidth * (FloorArray(2) - FloorArray(1)));
    end
end

% 3d plotting
% Start by convertin diff to rear ride height values
RRHarray = Diff2staticarray + 0.11;
% Creating an average ride height array
avgRideHeight = avgdist2road(1, :)';
% Create the meshgrid arrays for X and Y values
[avgRideHeightmesh, vCararraymesh] = meshgrid(avgRideHeight, vCararray);
figure
colormap(turbo)
scatter3(avgRideHeightmesh, vCararraymesh, DownforceVals, '.')
ylabel('vCar (kph)')
xlabel('Average ride height (m)')
zlabel('Downforce generated (N)')

%% Creating the fitted function
% Need to start by constructing the input vectors for the fit function
num = length(avgRideHeightmesh(1, :)) * length(avgRideHeightmesh(:, 1));
avgRideHeightcolumn = reshape(avgRideHeightmesh, [num, 1]);
vCarcolumn = reshape(vCararraymesh, [num, 1]);

% And the downforce vector
DownforceValscolumn = reshape(DownforceVals, [num, 1]);

% Fitting
fittedDWFfunc = fit([avgRideHeightcolumn vCarcolumn], DownforceValscolumn, 'poly55', 'normalize', 'on');

% Prediciting values, start by normalizing inputs
avgRideHeightcolumnnormal = (avgRideHeightcolumn - mean(avgRideHeightcolumn)) / std(avgRideHeightcolumn);
vCarcolumnnormal = (vCarcolumn - mean(vCarcolumn)) / std(vCarcolumn);
x = avgRideHeightcolumnnormal;
y = vCarcolumnnormal;

% Predicting
predictedvals = fittedDWFfunc.p00 + fittedDWFfunc.p10.*x + fittedDWFfunc.p01.*y + fittedDWFfunc.p20.*x.^2 + fittedDWFfunc.p11.*x.*y + fittedDWFfunc.p02.*y.^2 ... 
    + fittedDWFfunc.p30.*x.^3 + fittedDWFfunc.p21.*x.^2.*y + fittedDWFfunc.p12.*x.*y.^2 + fittedDWFfunc.p03.*y.^3 + fittedDWFfunc.p40.*x.^4 ...
    + fittedDWFfunc.p31.*x.^3.*y + fittedDWFfunc.p22.*x.^2.*y.^2 + fittedDWFfunc.p13.*x.*y.^3 + fittedDWFfunc.p04.*y.^4 + fittedDWFfunc.p50.*x.^5 ...
    + fittedDWFfunc.p41.*x.^4.*y + fittedDWFfunc.p32.*x.^3.*y.^2 + fittedDWFfunc.p23.*x.^2.*y.^3 + fittedDWFfunc.p14.*x.*y.^4 + fittedDWFfunc.p05.*y.^5;

% predictedvals = fittedDWFfunc([avgRideHeightcolumn vCarcolumn]);

predictedvalsmesh = reshape(predictedvals, [length(avgRideHeightmesh(:, 1)) length(avgRideHeightmesh(1, :))]);

% Plotting with scattered actual values
hold on
surf(avgRideHeightmesh, vCararraymesh, predictedvalsmesh)
