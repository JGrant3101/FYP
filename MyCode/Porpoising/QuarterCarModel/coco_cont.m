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
Inputs(4, 1) = 3400;
% Tyre vertical stiffness (N/m)
Inputs(5, 1) = 0.9*2.7 * 10^5;
% Static ride height (m)
Inputs(6, 1) = 0.1;
% vCar (kph)
Inputs(7, 1) = 280;
% Upper downforce elements multiplier
Inputs(8, 1) = 0.365;
% Mean for Inverse Gaussian distribution
Inputs(9, 1) = 0.0001;
% Shape factor for Inverse Gaussian distribution
Inputs(10, 1) = 2.4;
% Scaling applied to Inverse Gaussian distribution
Inputs(11, 1) = 0.31*(500/9)^2;

% Defining our initial x array
% x0 = [-31.1e-3;-8.4e-3;0;0];
x0 = [-0.0371;-0.01;0;0];
% x0 = [0;0;0;0];

% Defining the names of the inputs already written above
pnames = {'Ms' 'Mu' 'Ks' 'Cs' 'Kt' 'H' 'vCar' 'A' 'mew' 'lamda' 'scaling'};

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
Index = 7;
bdread0 = varyingparameters(param, [50 350], prob, args, Inputs);

%% Further inspection of these HB points

labs = coco_bd_labs(bdread0, 'HB');

HBbd = cell(1, length(labs));
Po = cell(1, length(labs));

for i = 1:length(labs)
    param4HB = 'H';
    Index4HB = 6;
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
    HBbd{i} = coco(prob1, sprintf('Test%d', i), [], {param, param4HB}, {[50 350], [0.05 0.12]});

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

    % Initialise figure for plotting of Zs and Zu evolutions
    fig2 = figure(i+2); ax2 = axes('Parent', fig2); hold on; view(3); xlabel(param); ylabel('Zs (m)'); zlabel(param4HB);title('Evolution of Zs with '+convertCharsToStrings(param)+' passing through a HB point') 
    % Plot Equilibrium results for Zs
    paramep = coco_bd_col(HBbd{1,1}, param); % Get the ep param at every HB point
    param4HBvals = coco_bd_col(HBbd{1,1}, param4HB); % Get the HB parameter at every HB point
    xHB = coco_bd_col(HBbd{1,1}, 'x'); % Get the x state at every HB point
    % lambdaep = coco_bd_col(bdread0, 'eigs'); % Get the x state
    % stabep = any(real(lambdaep)>0) ;
    % Plot Zs vs param varied in ep contin and param varied in HB contin
    plot3(ax2, paramep, xHB(1,:), param4HBvals, 'r');
    % plot3(ax2, paramep(stabep), xep(1,stabep), param4HBep(stabep), '.r');

    fig3 = figure(i+3); ax3 = axes('Parent', fig3); hold on; view(3); xlabel(param); ylabel('Zu (m)'); zlabel(param4HB);title('Evolution of Zu with '+convertCharsToStrings(param)+' passing through a HB point')
    % Plot Zu vs param varied in ep contin and param varied in HB contin
    plot3(ax3, paramep, xHB(2,:), param4HBvals, 'r');
    % plot3(ax3, paramep(stabep), xep(2,stabep), param4HBep(stabep), '.r') ;
    
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

            % Plot PO results for Zs
            bd2 = coco_bd_read(sprintf('po_run%d', i));
            parampo = coco_bd_col(bd2, param); % Get the current parameter
            param4HBpo = coco_bd_col(bd2, param4HB);
            x_max = coco_bd_col(bd2, 'x_max')';   % Get the state vector
            x_min = coco_bd_col(bd2, 'x_min')';   % Get the state vector
            stabpo = coco_bd_col(bd2, 'po.test.USTAB') == 0; % Get the stability
        
            plot3(ax2, parampo(stabpo), x_max(1,stabpo), param4HBpo(stabpo),'.k') ;
            plot3(ax2, parampo(~stabpo), x_max(1,~stabpo), param4HBpo(~stabpo),'.m') ;
            plot3(ax2, parampo(stabpo), x_min(1,stabpo), param4HBpo(stabpo),'.k') ;
            plot3(ax2, parampo(~stabpo), x_min(1,~stabpo), param4HBpo(~stabpo),'.m') ;
            % if any(~stabpo)
            %     legend('Stable ep', 'Unstable ep', 'Stable po', 'Unstable po')
            % else 
            %     legend('Stable ep', 'Unstable ep', 'Stable po')
            % end
    
            % Plot PO results for Zu
            plot3(ax3, parampo(stabpo), x_max(2,stabpo), param4HBpo(stabpo),'.k') ;
            plot3(ax3, parampo(~stabpo), x_max(2,~stabpo), param4HBpo(~stabpo),'.m') ;
            plot3(ax3, parampo(stabpo), x_min(2,stabpo), param4HBpo(stabpo),'.k') ;
            plot3(ax3, parampo(~stabpo), x_min(2,~stabpo), param4HBpo(~stabpo),'.m') ;
            % if any(~stabpo)
            %     legend('Stable ep', 'Unstable ep', 'Stable po', 'Unstable po')
            % else 
            %     legend('Stable ep', 'Unstable ep', 'Stable po')
            % end
        end
    end
   legend(ax2, {'Stable ep', 'Unstable ep', 'Stable po', 'Unstable po'})
   legend(ax3, {'Stable ep', 'Unstable ep', 'Stable po', 'Unstable po'})
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

