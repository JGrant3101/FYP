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
Inputs(7, 1) = 200;
% Upper downforce elements multiplier
Inputs(8, 1) = 0.365;
% Mean for Inverse Gaussian distribution
Inputs(9, 1) = 0.033;
% Shape factor for Inverse Gaussian distribution
Inputs(10, 1) = 0.07;
% Scaling applied to Inverse Gaussian distribution
Inputs(11, 1) = 0.1*(500/9)^2;

u0 = [-31.1e-3;-8.4e-3;0;0] ;
prob = coco_prob();

bd1 = coco(prob, 'car/4', 'ode', 'isol', 'ep', ...
  @(u,p)Suspension(0,u,p), [], [], u0, {'al' 'be' 'de' 'ro'}, Inputs,  ...               % ep toolbox arguments
  1, {'be' 'ep.test.SN' 'ep.test.HB' 'ep.test.USTAB' 'atlas.test.FP'}, [2 7]); % cont toolbox arguments