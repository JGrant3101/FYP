%% Testing the quarter car model by running time series data through it
% Start by defining the time array and corresponding velocity array
% Want to accelerate from 0 to 320kph at 5ms^-2 to start
vCar = [0:0.5:320]';
t = [0:0.1:(length(vCar)/10)-0.1]';
% Setting up the constant input terms for the suspension function
Inputs = [];
% Sprung mass (kg)
Inputs(1, 1) = 180;
% Unsprung mass (kg)
Inputs(2, 1) = 50;
% Suspension stiffness (N/m)
Inputs(3, 1) = 3*10^4;
% Suspension damping (Ns/m)
Inputs(4, 1) = 2000;
% Tyre vertical stiffness (N/m)
Inputs(5, 1) = 2.7 * 10^5;
% Static ride height (mm)
Inputs(6, 1) = 100;
% vCar, this will be populated later
Inputs(7, 1) = 200;
% Upper downforce elements multiplier
Inputs(8, 1) = 0.365;
% Mean for Inverse Gaussian distribution
Inputs(9, 1) = 35;
% Shape factor for Inverse Gaussian distribution
Inputs(10, 1) = 8;
% Scaling applied to Inverse Gaussian distribution
Inputs(11, 1) = 2.3*(500/9)^2;


sol = ode45(@(t, x)Suspension(t, x, Inputs), [0, 64], [-0.031; -0.0084; 0; 0], odeset('RelTol', 1e-8)); % Simulate 


return