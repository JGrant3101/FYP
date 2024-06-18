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
Inputs(3, 1) = 0.9*10^5;
% Suspension damping (Ns/m)
Inputs(4, 1) = 3775;
% Tyre vertical stiffness (N/m)
Inputs(5, 1) = 0.9*2.7 * 10^5;
% Static ride height (m)
Inputs(6, 1) = 0.1;
% vCar (kph)
Inputs(7, 1) = 280;
% Upper downforce elements multiplier
Inputs(8, 1) = 0.365;
% Mean for Log Normal Distribution
Inputs(9, 1) = 0.0001;
% Shape factor for Log Normal Distribution
Inputs(10, 1) = 2.4;
% Scaling applied to Log Normal Distribution
Inputs(11, 1) = 0.31*(500/9)^2;

% Defining our initial x array
x0 = [-0.0371;-0.01;0;0];

% Defining the names of the inputs already written above
pnames = {'Ms' 'Mu' 'Ks' 'Cs' 'Kt' 'H' 'vCar' 'A' 'mew' 'lamda' 'scaling'};

% Initialising the COCO problem
prob = coco_prob();
prob = coco_set(prob, 'ep', 'NSA', true);

% Defining the functions
SystemSetup = {@Suspension};

% Defining the function arguments
args = {SystemSetup{:}, x0, pnames, Inputs};

% Change the maximum step size
prob = coco_set(prob, 'cont', 'h_max', 1);

prob = coco_set(prob, 'cont', 'h_min', 1e-10);

% Increase the number of continuation steps
prob = coco_set(prob, 'cont', 'PtMX', 500);
% Increase the number of Iterations
%prob = coco_set(prob, 'cont', 'ItMX', 50);


%% Varying single parameter to find HB points

% Calling function to vary a parameter and run a continuation
param = 'vCar';
Index = 7;
bdread0 = varyingparameters(param, [50 350], prob, args, Inputs);

%% Plotting the HB point of the initial system
figure; clf; hold on
thm1 = struct('special', {{'HB', 'EP'}});
coco_plot_bd(thm1, 'Initial', 'vCar', '||x||_2')
title('Graph of the equilibrium point of the system as car speed varies')
fontsize(gca, 20, 'points')
grid on
%% Further inspection of these HB points

labs = coco_bd_labs(bdread0, 'HB');

HBbd = cell(1, length(labs));
Po = cell(1, length(labs));

for i = 1:length(labs)
    param4HB = 'Ks';
    Index4HB = 3;
    % Creating an empty COCO problem structure
    prob1 = coco_prob();
    % Defining the settings of the problem structure
    prob1 = coco_set(prob1, 'ep', 'NSA', true);
    prob1 = coco_set(prob1, 'cont', 'h_max', 1000);
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
    arrayparamep = linspace(minparamep, maxparamep, 11);

    % Adding an event to get points of interest
    prob1 = coco_add_event(prob1, 'POI', param, arrayparamep);

    % Running COCO
    HBbd{i} = coco(prob1, sprintf('Test%d', i), [], {param, param4HB}, {[50 350], [10000 400000]});

    % Plotting the HB2HB continuation results
    figure(i+1); clf; hold on
    thm1 = struct('special', {{'HB', 'EP'}});
    thm2 = struct('special', {{'EP'}});
    % coco_plot_bd(thm1, 'Initial', param, param4HB, '||x||_2')
    coco_plot_bd(thm2, sprintf('Test%d', i), param4HB, param)
    ylabel('vCar (kph)')
    xlabel('Ks (N/m)')

    % Want to also scatter the original point tested to get a reference of
    % where that is in relation to our HB points
    row = find(cell2mat([bdread0(2:end,1)]) == 0);
    row = row(1) + 1;
    [~,column] = find(strcmp(bdread0,'||x||_2'));
    scatter(Inputs(Index4HB, 1), 300.73, 100, 'green', 'square', 'filled')
    title(['Family of HB points for varying ' param4HB])
    hold off; grid on; fontsize(gca, 20, 'points')
    
    % Want to plot a series of curves that show the values of Zs and Zu
    % around the equilibirum after the HB point
    % Start by finding the array of labs of points of interest
    HBlabs = coco_bd_labs(HBbd{1,1}, 'POI');
    
    % Run a for loop within this for loop to get the HB2PO continuations
    % for each of the chosen HB points
    for j = 1:length(HBlabs)
        prob2 = coco_prob();
        prob2 = coco_set(prob2, 'coll', 'NTST', 1000);
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

    x_val = cell2mat([bdread0(2:end,28)]') ;
    figure ; plot(IndependentVar_val, x_val(1:2,:)) ;
    hep = Inputs(6, 1) + x_val(1,:) + x_val(2,:);
    hold on
    plot(IndependentVar_val, hep)
    xlabel(param2vary)
    ylabel('Equilibrium values (m)')
    legend('Zs (m)', 'Zu (ms)', 'Ride height (m)')
    fontsize(gca, 20, 'points')
end

