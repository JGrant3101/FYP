clear all
close all
clc

%% Initialising inputs

% Setting up the constant input terms for the suspension function
Inputs = [];
% Unsprung front mass (kg)
Inputs(1, 1) = 50;
% Unsprung rear mass (kg)
Inputs(2, 1) = 50;
% Sprung mass (kg)
Inputs(3, 1) = 600;
% Sprung mass moment of inertia (kgm^2)
Inputs(4, 1) = 670;
% Distance from front axle to CoG (m)
Inputs(5, 1) = 1.98;
% Distance from rear axle to CoG (m)
Inputs(6, 1) = 1.62;
% Front tyre vertical stiffness (N/m)
Inputs(7, 1) = 2 * 2.7 * 10^5;
% Front Suspension stiffness (N/m)
Inputs(8, 1) = 1.5 * 2 * 10^5;
% Rear tyre vertical stiffness (N/m)
Inputs(9, 1) =  2 * 0.9 * 2.7 * 10^5;
% Rear Suspension stiffness (N/m)
Inputs(10, 1) = 1.5 * 1.81 * 10^5;
% Front damping (Ns/m)
Inputs(11, 1) = 9247;
% Rear damping (Ns/m)
Inputs(12, 1) = 9739;
% Static front ride height (m)
Inputs(13, 1) = 0.0943;
% Static rear ride height (m)
Inputs(14, 1) = 0.11;
% Car speed (kph)
Inputs(15, 1) = 270;
% Front wing aero multiplier
Inputs(16, 1) = 0.254;
% Rear wing aero multiplier
Inputs(17, 1) = 0.636;
% Front downforce from floor gradient
Inputs(18, 1) = 0.5;
% Front downforce from floor constant
Inputs(41, 1) = 0.3;

% Reading in the mat file that contains the parameters needed for the floor
% aero model
load("AeroValues.mat")

% Downforce from floor polynomial terms
Inputs(19, 1) = fittedDWFfunc.p00;
Inputs(20, 1) = fittedDWFfunc.p10;
Inputs(21, 1) = fittedDWFfunc.p01;
Inputs(22, 1) = fittedDWFfunc.p20;
Inputs(23, 1) = fittedDWFfunc.p11;
Inputs(24, 1) = fittedDWFfunc.p02;
Inputs(25, 1) = fittedDWFfunc.p30;
Inputs(26, 1) = fittedDWFfunc.p21;
Inputs(27, 1) = fittedDWFfunc.p12;
Inputs(28, 1) = fittedDWFfunc.p03;
%Inputs(29, 1) = fittedDWFfunc.p40;
Inputs(29, 1) = fittedDWFfunc.p31;
Inputs(30, 1) = fittedDWFfunc.p22;
Inputs(31, 1) = fittedDWFfunc.p13;
Inputs(32, 1) = fittedDWFfunc.p04;
%Inputs(34, 1) = fittedDWFfunc.p50;
%Inputs(35, 1) = fittedDWFfunc.p41;
Inputs(33, 1) = fittedDWFfunc.p32;
Inputs(34, 1) = fittedDWFfunc.p23;
Inputs(35, 1) = fittedDWFfunc.p14;
Inputs(36, 1) = fittedDWFfunc.p05;
Inputs(37, 1) = avgRHmean;
Inputs(38, 1) = avgRHstd;
Inputs(39, 1) = vCarmean;
Inputs(40, 1) = vCarstd;

% Floor road interaction parameters
% Coefficient of resitution
Inputs(42, 1) = 0.8;

% Defining our initial x array
% x0 = [-31.1e-3;-8.4e-3;0;0];
% x0 = [-0.01; -0.01; -0.02; 0.0044; 0; 0; 0; 0];
x0 = [-0.02; -0.03; -0.07; -0.009; 0; 0; 0; 0];
% x0 = [0;0;0;0];

% Defining the names of the inputs already written above
pnames = {'MuF' 'MuR' 'Ms' 'Is' 'LF' 'LR' 'KtF' 'KsF' 'KtR' 'KsR' 'CsF' 'CsR' 'StaticFRH' 'StaticRRH' 'vCar' 'AF' 'AR' 'm' 'p00' 'p10' ...
    'p01' 'p20' 'p11' 'p02' 'p30' 'p21' 'p12' 'p03' 'p31' 'p22' 'p13' 'p04' 'p32' 'p23' 'p14' 'p05' 'avgRHmean' 'avgRHstd' 'vCarmean' ...
    'vCarstd' 'DWFFloorFConstant' 'e'};

% Initialising the COCO problem
prob = coco_prob();
prob = coco_set(prob, 'ep', 'NSA', true);

% Defining the functions
SystemSetup = {@Suspension}; %, @Suspension_dx, @Suspension_dp};

% Defining the function arguments
args = {SystemSetup{:}, x0, pnames, Inputs};

% Change the maximum step size
prob = coco_set(prob, 'cont', 'h_max', 1);

prob = coco_set(prob, 'cont', 'h_min', 1e-10);

% Increase the number of continuation steps
prob = coco_set(prob, 'cont', 'PtMX', 500);
% Increase the number of Iterations
%prob = coco_set(prob, 'cont', 'ItMX', 50);

% Adding the lyapunov function 
% [data, uidx] = coco_get_func_data(prob, 'ep', 'data', 'uidx');
% 
% prob = coco_add_func(prob, 'lyap', @lyapunov, data.ep_eqn, 'regular', 'L1', 'uidx', uidx);

%% Varying single parameter to find HB points

% Calling function to vary a parameter and run a continuation
param = 'vCar';
Index = 15;
bdread0 = varyingparameters(param, [50 350], prob, args, Inputs);

%% Further inspection of these HB points

labs = coco_bd_labs(bdread0, 'HB');

HBbd = cell(1, length(labs));
Po = cell(1, length(labs));

for i = 1:length(labs)
    param4HB = 'A';
    Index4HB = 8;
    % Creating an empty COCO problem structure
    prob1 = coco_prob();
    % Defining the settings of the problem structure
    prob1 = coco_set(prob1, 'ep', 'NSA', true);
    prob1 = coco_set(prob1, 'cont', 'h_max', 1);
    prob1 = coco_set(prob1, 'cont', 'h_min', 1e-10);
    prob1 = coco_set(prob1, 'cont', 'PtMX', 2000);
    prob1 = coco_set(prob1, 'cont', 'ItMX', 5000);
    prob1 = coco_set(prob1, 'cont', 'NAdapt', 1); 

    % Restarting the continuation along a family of Hopf Bifurcations
    prob1 = ode_HB2HB(prob1, '', 'Initial', labs(i));

    % Adding the lyapunov function 
    [data1, uidx1] = coco_get_func_data(prob1, 'ep', 'data', 'uidx');
    
    prob1 = coco_add_func(prob1, 'lyap', @lyapunov, data1.ep_eqn, 'regular', 'L1', 'uidx', uidx1);

    % Adding an event to check where the lyapunov coefficient is 0
    prob1 = coco_add_event(prob1, 'GH', 'L1', 0);

    % Creating the custom POI event to identify HB points that in the
    % future will have a PO contin run on them
    minparamep = min(coco_bd_col(bdread0, param));
    maxparamep = max(coco_bd_col(bdread0, param));
    % minparamep = 286;
    % maxparamep = 310;
    arrayparamep = linspace(minparamep, maxparamep, 11);

    % Adding an event to get points of interest
    prob1 = coco_add_event(prob1, 'POI', param, arrayparamep);

    % prob1 = coco_add_slot(prob1, 'slot_bd_min_max', @slot_bd_min_max, [], 'bddat');

    % Running COCO
    HBbd{i} = coco(prob1, sprintf('Test%d', i), [], {param, param4HB}, {[50 350], [0.2 2]});

    % Plotting the HB2HB continuation results
    figure(i+1); clf; hold on
    thm1 = struct('special', {{'HB', 'EP'}});
    thm2 = struct('special', {{'EP'}});
    coco_plot_bd(thm1, 'Initial', param, param4HB, '||x||_2')
    coco_plot_bd(thm2, sprintf('Test%d', i), param, param4HB, '||x||_2')

    % Want to also scatter the original point tested to get a reference of
    % where that is in relation to our HB points
    row = find(cell2mat([bdread0(2:end,1)]) == 0);
    row = row(1) + 1;
    [~,column] = find(strcmp(bdread0,'||x||_2'));
    scatter3(Inputs(Index, 1), Inputs(Index4HB, 1), bdread0{row, column}, 100, 'green', 'square', 'filled')
    hold off; grid on; view(3)
    
    % Want to plot a series of curves that show the values of Zs and Zu
    % around the equilibirum after the HB point
    % Start by finding the array of labs of points of interest
    HBlabs = coco_bd_labs(HBbd{1,1}, 'POI');
    
    % Run a for loop within this for loop to get the HB2PO continuations
    % for each of the chosen HB points
    for j = 1:length(HBlabs)
        prob2 = coco_prob();
        prob2 = coco_set(prob2, 'coll', 'NTST', 500);
        prob2 = coco_set(prob2, 'po', 'bifus', true);
        prob2 = ode_HB2po(prob2, '', sprintf('Test%d', i), HBlabs(j));
        prob2 = coco_set(prob2, 'cont', 'PtMX', 200);
        prob2 = coco_set(prob2, 'cont', 'ItMX', 1500);
        prob2 = coco_set(prob2, 'cont', 'h_min', 1e-10);
        prob2 = coco_set(prob2, 'cont', 'NAdapt', 1, 'h_max', 1);
        % 
        % % Tell COCO to store extra information about the periodic orbits
        prob2 = coco_add_slot(prob2, 'slot_bd_min_max', @slot_bd_min_max, [], 'bddat');
    
        Po{i} = coco(prob2, sprintf('po_run%d', i), [], 1, {param 'po.period'}, [50 350]);
        
        if Po{1, 1}{2, 1} ~= -1
            % Updating inputs before running an EP about the new point for
            % plotting
            Inputsep = Inputs;
            Inputsep(7, 1) = 20;
            Inputsep(Index4HB, 1) = Po{1, 1}{2, Index4HB + 8};

            args = {SystemSetup{:}, x0, pnames, Inputsep};

            variableargs = {1, 'vCar', [10 350]};

            % Running a relevant EP for the parameter that was changed
            bdep = coco(prob, sprintf('epfrompo_run%d', j), @ode_isol2ep, args{:}, variableargs{:});
            
            % Plotting Zs results
            % Plotting the ep values
            figure
            xep = coco_bd_col(bdep, 'x')';
            zsep = xep(:, 1);
            vCarep = coco_bd_col(bdep, 'vCar')';

            plot(vCarep, zsep)
            hold on

            % Plotting the PO values
            vCarPO = coco_bd_col(Po{i}, param);
            x_max = coco_bd_col(Po{i}, 'x_max')';   % Get the state vector
            x_min = coco_bd_col(Po{i}, 'x_min')';   % Get the state vector
            stabpo = coco_bd_col(Po{i}, 'po.test.USTAB') == 0; % Get the stability
            
            zs_max = x_max(1, :);
            zs_min = x_min(1, :);

            plot(vCarPO(stabpo), zs_max(stabpo),'.k') ;
            plot(vCarPO(~stabpo), zs_max(~stabpo),'.m') ;
            plot(vCarPO(stabpo), zs_min(stabpo),'.k') ;
            plot(vCarPO(~stabpo), zs_min(~stabpo),'.m') ;

            xlabel('vCar (kph)')
            ylabel('Zs (m)')
            title(['EP and PO values of Zs across vCar range with ' param4HB ' set to: ' num2str(Inputsep(Index4HB, 1))])
            
            legend('Ep values', 'Stable PO results', 'Unstable PO results')
            % Plotting the Zu results
            % Plotting the ep values
            figure
            zuep = xep(:, 2);

            plot(vCarep, zuep)
            hold on

            % Plotting the PO values
            zu_max = x_max(2, :);
            zu_min = x_min(2, :);

            plot(vCarPO(stabpo), zu_max(stabpo),'.k') ;
            plot(vCarPO(~stabpo), zu_max(~stabpo),'.m') ;
            plot(vCarPO(stabpo), zu_min(stabpo),'.k') ;
            plot(vCarPO(~stabpo), zu_min(~stabpo),'.m') ;

            xlabel('vCar (kph)')
            ylabel('Zu (m)')
            title(['EP and PO values of Zu across vCar range with ' param4HB ' set to: ' num2str(Inputsep(Index4HB, 1))])

            legend('Ep values', 'Stable PO results', 'Unstable PO results')
        end
    end
end


%% Section for functions
function bdread0 = varyingparameters(param2vary, range, prob, args, Inputs)
    variableargs = {1, param2vary, range};

    bd0 = coco(prob, 'Initial', @ode_isol2ep, args{:}, variableargs{:});
    
    
    bdread0 = coco_bd_read('Initial'); % Get the Omega parameter
    [~,column] = find(strcmp(bdread0,param2vary));
    IndependentVar_val = cell2mat([bdread0(2:end,column)])' ;

    x_val = cell2mat([bdread0(2:end,59)]') ;
    figure ; plot(IndependentVar_val, x_val(1:2,:)) ;
    hep = Inputs(6, 1) + x_val(1,:) + x_val(2,:);
    hold on
    plot(IndependentVar_val, hep)
    xlabel(param2vary)
    ylabel('Equilibrium values (m)')
    legend('Zs (m)', 'Zu (ms)', 'Ride height (m)')
end

