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

%% Further inspection of these HB points

labs = coco_bd_labs(bdread0, 'HB');

HBbd = cell(1, length(labs));
Po = cell(1, length(labs));

for i = 1:length(labs)
    param4HB = 'scaling';
    Index4HB = 11;
    % Creating an empty COCO problem structure
    prob1 = coco_prob();
    % Defining the settings of the problem structure
    prob1 = coco_set(prob1, 'ep', 'NSA', true);
    prob1 = coco_set(prob1, 'cont', 'h_max', 2);
    prob1 = coco_set(prob1, 'cont', 'h_min', 1e-10);
    prob1 = coco_set(prob1, 'cont', 'PtMX', 700);
    prob1 = coco_set(prob1, 'cont', 'ItMX', 500);
    prob1 = coco_set(prob1, 'cont', 'NAdapt', 4); 

    % Restarting the continuation along a family of Hopf Bifurcations
    prob1 = ode_HB2HB(prob1, '', 'Initial', labs(i));

   % prob1 = coco_add_slot(prob1, 'slot_bd_min_max', @slot_bd_min_max, [], 'bddat');

    % Running COCO
    HBbd{i} = coco(prob1, sprintf('Test%d', i), [], {param, param4HB}, {[50 350], [400 1600]});
    
    prob2 = coco_prob();
    prob2 = coco_set(prob2, 'coll', 'NTST', 200);
    prob2 = coco_set(prob2, 'po', 'bifus', true);
    prob2 = ode_HB2po(prob2, '', 'Initial', labs(i));
    prob2 = coco_set(prob2, 'cont', 'PtMX', 1500);
    prob2 = coco_set(prob2, 'cont', 'ItMX', 500);
    prob2 = coco_set(prob2, 'cont', 'h_min', 1e-10);
    prob2 = coco_set(prob2, 'cont', 'NAdapt', 1, 'h_max', 1);

    % Tell COCO to store extra information about the periodic orbits
    prob2 = coco_add_slot(prob2, 'slot_bd_min_max', @slot_bd_min_max, [], 'bddat');

    Po{i} = coco(prob2, sprintf('po_run%d', i), [], 1, {param 'po.period'}, [50 350]);

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

    figure(i+2); clf; hold on; xlabel(param); ylabel('Zs (m)'); title('Evolution of Zs with ' + convertCharsToStrings(param) + ' passing through a HB point') 
    % Plot Equilibrium results for Zs
    paramep = coco_bd_col(bdread0, param); % Get the f parameter
    xep = coco_bd_col(bdread0, 'x'); % Get the x state
    lambdaep = coco_bd_col(bdread0, 'eigs'); % Get the x state
    stabep = any(real(lambdaep)>0) ;
    plot(paramep(~stabep), xep(1,~stabep), '.b') ;
    plot(paramep(stabep), xep(1,stabep), '.r') ;

    % Plot PO results for Zs
    bd2 = coco_bd_read(sprintf('po_run%d', i));
    parampo = coco_bd_col(bd2, param); % Get the current parameter
    x_max = coco_bd_col(bd2, 'x_max')';   % Get the state vector
    x_min = coco_bd_col(bd2, 'x_min')';   % Get the state vector
    stabpo = coco_bd_col(bd2, 'po.test.USTAB') == 0; % Get the stability

    plot(parampo(stabpo), x_max(1,stabpo),'.k') ;
    plot(parampo(~stabpo), x_max(1,~stabpo),'.m') ;
    plot(parampo(stabpo), x_min(1,stabpo),'.k') ;
    plot(parampo(~stabpo), x_min(1,~stabpo),'.m') ;
    if any(~stabpo)
        legend('Stable ep', 'Unstable ep', 'Stable po', 'Unstable po')
    else 
        legend('Stable ep', 'Unstable ep', 'Stable po')
    end

    figure(i+3); clf; hold on; xlabel(param); ylabel('Zu (m)'); title('Evolution of Zu with ' + convertCharsToStrings(param) + ' passing through a HB point')
    % Plot Equilibrium results for Zu
    plot(paramep(~stabep), xep(2,~stabep), '.b') ;
    plot(paramep(stabep), xep(2,stabep), '.r') ;

    % Plot PO results for Zu
    plot(parampo(stabpo), x_max(2,stabpo),'.k') ;
    plot(parampo(~stabpo), x_max(2,~stabpo),'.m') ;
    plot(parampo(stabpo), x_min(2,stabpo),'.k') ;
    plot(parampo(~stabpo), x_min(2,~stabpo),'.m') ;
    if any(~stabpo)
        legend('Stable ep', 'Unstable ep', 'Stable po', 'Unstable po')
    else 
        legend('Stable ep', 'Unstable ep', 'Stable po')
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
end

