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
TunnelWidthFraction = 0.3;
TunnelWidth = FloorWidth * TunnelWidthFraction;
CondiWidth = (1 - TunnelWidthFraction) * FloorWidth;
FloorIndices = 30001;

% Create an array of floor distances based on the floor length
FloorArray = linspace(0, FloorLength, FloorIndices)';

% Start by defining static front and rear ride heights for the condi
% section
StaticRRH = 0.15;
FloorAngle = 0.25;
StaticFRH = StaticRRH - FloorLength * sind(FloorAngle);
%StaticFRH = 0.0943;


% Defining the x and y for the low point of the floor's condi section
Xlow = 1.7;
Ylow = 0.12;
% Will form a quadratic defining the floor profile
temp1 = zeros([4, 4]);
temp2 = zeros([4, 1]);

% Lowering the floor here while working on code
StaticFRH = StaticFRH;
StaticRRH = StaticRRH;
Ylow = Ylow;

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

% Now defining the tunnel section
StaticRearMin = 0.08;
StaticFrontMin = StaticRearMin - FloorLength * sind(FloorAngle);

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
vCar = 120; % This is in kph will convert to m/s
vCar = (vCar * 10^3)/3600;

% Calculating inlet mass flow
mdot = density * vCar * FloorWidth * StaticFRH;

% Splitting to the two sections
mdotcondi = 0.9 * mdot;
mdotTunnel = 0.1 * mdot;

% Recalculating airspeed for the condi and tunnel sections
airspeedcondi = mdotcondi / (CondiWidth * StaticFRH * density);
airpseedTunnel = mdotTunnel / (TunnelWidth * StaticFrontMin * density);

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