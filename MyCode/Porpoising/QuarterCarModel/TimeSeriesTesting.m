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
Inputs(4, 1) = 200;
% Tyre vertical stiffness (N/m)
Inputs(5, 1) = 2.7 * 10^5;
% Static ride height (mm)
Inputs(6, 1) = 100;
% vCar, this will be populated later
Inputs(7, 1) = 300;
% Upper downforce elements multiplier
Inputs(8, 1) = 0.365;
% Mean for Inverse Gaussian distribution
Inputs(9, 1) = 35;
% Shape factor for Inverse Gaussian distribution
Inputs(10, 1) = 8;
% Scaling applied to Inverse Gaussian distribution
Inputs(11, 1) = 2.3*(500/9)^2;


%sol = ode45(@(t, x)Suspension(t, x, Inputs), [0, 100], 0.01*rand(4, 1), odeset('RelTol', 1e-8)); % Simulate


%return
% Now that all inputs have been set up can run the function in a for loop
% (which is required as ride heights calculated at the previous time step will have an impact on those at the next)
% to calculate the ride heights are each speed value.
% First want to initialise an x vector of all zeros
x = zeros(4, 1);
% Then want to initialise arrays which will be populated with results
Zs = zeros(length(t), 1);
Zu = zeros(length(t), 1);
RideHeight = zeros(length(t), 1);
RideHeight(:, 1) = 120;
DWFUpper = zeros(length(t), 1);
DWFFloor = zeros(length(t), 1);
for i = 2:length(t)
    % Assinging vCar
    Inputs(7, 1) = vCar(i, 1);
    % Calling the suspension function
    [rhs, VarOfInterest] = Suspension(x, Inputs);
    % Saving variables of interest
    Zs(i, 1) = Zs(i-1, 1) + 0.1 * rhs(1, 1); 
    Zu(i, 1) = Zu(i-1, 1) + 0.1 * rhs(2, 1);
    DWFUpper(i, 1) = VarOfInterest(1, 1);
    DWFFloor(i, 1) = VarOfInterest(2, 1);
    
    % Creating new x array
    x(1, 1) = Zs(i, 1);
    x(2, 1) = Zu(i, 1);
    x(3, 1) = rhs(1, 1);
    x(4, 1) = rhs(2, 1);
end

RideHeight = RideHeight - Zs - Zu;
