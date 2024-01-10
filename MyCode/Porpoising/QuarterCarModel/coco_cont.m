clear all
close all
clc

% Setting up the constant input terms for the suspension function
Inputs = [];
% Sprung mass (kg)
Inputs(1, 1) = 180;
% Unsprung mass (kg)
Inputs(2, 1) = 50;
% Suspension stiffness (N/m)
Inputs(3, 1) = 10^5;
% Suspension damping (Ns/m)
Inputs(4, 1) = 3400;
% Tyre vertical stiffness (N/m)
Inputs(5, 1) = 2.7 * 10^5;
% Static ride height (m)
Inputs(6, 1) = 0.1;
% vCar (kph)
Inputs(7, 1) = 250;
% Upper downforce elements multiplier
Inputs(8, 1) = 0.365;
% Mean for Inverse Gaussian distribution
Inputs(9, 1) = 0.033;
% Shape factor for Inverse Gaussian distribution
Inputs(10, 1) = 0.07;
% Scaling applied to Inverse Gaussian distribution
Inputs(11, 1) = 0.31*(500/9)^2;

% Defining our initial x array
x0 = [-31.1e-3;-8.4e-3;0;0];

% Defining the names of the inputs already written above
pnames = {'Ms' 'Mu' 'Ks' 'Cs' 'Kt' 'H' 'vCar' 'A' 'mew' 'lamda' 'scaling'};

% Initialising the COCO problem
prob = coco_prob();
prob = coco_set(prob, 'ep', 'NSA', true);

% Defining the functions
SystemSetup = {@Suspension, @Suspension_dx, @Suspension_dp};

% Defining the function arguments
args = {SystemSetup{:}, x0, pnames, Inputs};

% Change the maximum step size
prob = coco_set(prob, 'cont', 'h_max', 1);

prob = coco_set(prob, 'cont', 'h_min', 1e-10);

% Increase the number of continuation steps
prob = coco_set(prob, 'cont', 'PtMX', 300);
% Increase the number of Iterations
%prob = coco_set(prob, 'cont', 'ItMX', 50);

% Defining the arguments we want COCO to vary and in what range we want
% them varied
% variableargs = {2, {'Cs', 'vCar'}, {[2000 5000], [100 340]}};
variableargs = {1, 'vCar', [100 340]};

bd0 = coco(prob, 'Intitial', @ode_isol2ep, args{:}, variableargs{:});

bd = coco_bd_read('Intitial');

labs = coco_bd_labs(bd, 'HB');

prob1 = coco_prob();

prob1 = coco_set(prob1, 'cont1', 'PtMX', 50);

prob1 = ode_HB2HB(prob, '', 'Intitial', labs(2));

% Running COCO
bd1 = coco(prob, 'Test1', @ode_isol2ep, args{:}, variableargs{:});


%bd1 = coco(prob, 'car', 'ode', 'isol', 'ep', ...
 % @(u,p)Suspension(u,p), [], [], u0, {'al' 'be' 'de' 'ro'}, Inputs,  ...               % ep toolbox arguments
  %1, {'be' 'ep.test.SN' 'ep.test.HB' 'ep.test.USTAB' 'atlas.test.FP'}, [2 7]); % cont toolbox arguments






