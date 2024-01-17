clear all
close all
clc

%% Initialising inputs

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
Inputs(7, 1) = 300;
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
prob = coco_set(prob, 'cont', 'h_max', 200);

prob = coco_set(prob, 'cont', 'h_min', 1e-10);

% Increase the number of continuation steps
prob = coco_set(prob, 'cont', 'PtMX', 600);
% Increase the number of Iterations
%prob = coco_set(prob, 'cont', 'ItMX', 50);

%% Varying single parameter to find HB points

% Calling function to vary a parameter and run a continuation
bdread0 = varyingparameters('Kt', [100000 400000], prob, args, Inputs);

%% Further inspection of these HB points

labs = coco_bd_labs(bdread0, 'HB');

HBbd = cell(1, length(labs));

for i = 1:length(labs)
    
    % Creating an empty COCO problem structure
    prob1 = coco_prob();
    % Defining the settings of the problem structure
    prob1 = coco_set(prob1, 'ep', 'NSA', true);
    prob1 = coco_set(prob1, 'cont', 'h_max', 800);
    prob1 = coco_set(prob1, 'cont', 'h_min', 1e-10);
    prob1 = coco_set(prob1, 'cont', 'PtMX', 600);
    prob1 = coco_set(prob1, 'cont', 'ItMX', 500);
    prob1 = coco_set(prob1, 'cont', 'NAdapt', 4); 

    % Restarting the continuation along a family of Hopf Bifurcations
    prob1 = ode_HB2HB(prob1, '', 'Initial', labs(i));

    % Running COCO
    HBbd{i} = coco(prob1, sprintf('Test%d', i), [], {'Kt', 'Ks'}, {[100000 400000], [33000 220000]});

    figure(i+1); clf; hold on
    thm1 = struct('special', {{'HB', 'EP'}});
    thm2 = struct('sepcial', {{'EP'}});
    coco_plot_bd(thm1, 'Initial', 'Kt', 'Ks', '||x||_2')
    coco_plot_bd(thm2, sprintf('Test%d', i), 'Kt', 'Ks', '||x||_2')
    hold off; grid on; view(3)

end


%% Section for functions
function bdread0 = varyingparameters(param2vary, range, prob, args, Inputs)
    variableargs = {1, param2vary, range};

    bd0 = coco(prob, 'Initial', @ode_isol2ep, args{:}, variableargs{:});
    
    
    bdread0 = coco_bd_read('Initial'); % Get the Omega parameter
    [~,column] = find(strcmp(bdread0,param2vary));
    IndependentVar_val = cell2mat([bdread0(2:end,column)])' ;

    x_val = cell2mat([bdread0(2:end,28)]') ;
    figure ; plot(IndependentVar_val, x_val(1:2,:)) ;
    hep = Inputs(6, 1) + x_val(1,:) + x_val(2,:);
    hold on
    plot(IndependentVar_val, hep)
    xlabel(param2vary)
    ylabel('Equilibrium values (m)')
    legend('Zs (m)', 'Zu (ms)', 'Ride height (m)')
end

