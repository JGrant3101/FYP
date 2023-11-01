% DUFFING_DEMO Run a series of continuations on a Duffing equation

% Written by David A.W. Barton (david.barton@bristol.ac.uk) 2015
% Adapted by Ludovic Renson (l.renson@imperial.ac.uk) 2023

% Create a figure for plotting
figure;    
ax = gca();
hold(ax, 'on');
xlabel(ax, '\Omega');
ylabel(ax, '\Gamma');
zlabel(ax, '||x||');
view(ax, 3);

%% Generate an initial solution by simulation

p0 = [0.01, 1, 0.1, 0.1, 0.4]'; % Initial parameter values
%sol = ode45(@(t, x)duffing(t, x, p0), [0, 10], zeros(2, 1), odeset('RelTol', 1e-8)); % Simulate
t0 = linspace(0, 1, 101)'; % 101 data points initially
%x0 = deval(sol, sol.x(end) - t0(end) + t0)'; % Interpolate to get the required data

x0 = zeros(length(t0), 2);

%% Do an initial continuation in p2 to generate starting points for further continuation

% Create an empty COCO problem structure
prob = coco_prob();

% Tell COCO that it is a non-autonomous (time varying) problem
prob = coco_set(prob, 'ode', 'autonomous', false);

% Increase the number of mesh points used
prob = coco_set(prob, 'coll', 'NTST', 60);

% Turn adaptive meshing on (adapt every n steps)
prob = coco_set(prob, 'cont', 'NAdapt', 5); % Adapt every 5 steps

% Decrease the maximum step size
prob = coco_set(prob, 'cont', 'h_max', 1);

% Create a continuation problem: start with an isol (initial solution -
% user provided) and do a po continuation (periodic orbit), hence use
% the "isol2po" function
prob = ode_isol2po(prob, ...
    '', ...           % OID - only used if embedding this problem in a larger problem
    @duffing, ...     % Right-hand side of equations of motion (optionally the derivatives can also be provided to improve speed)
    t0, ...           % Mesh for the initial solution
    x0, ...           % Initial solution
    {'xi', 'omega', 'beta', 'Gamma', 'Omega'}, ... % Parameter names
    p0);              % Initial parameter values

% Add a user defined event to output the values at certain parameter values
prob = coco_add_event(prob, 'UZ', 'Gamma', [0.05, 0.1, 0.15, 0.2]); % Output when Gamma = 0.05, 0.1 or 0.2

% Tell COCO to store extra information about the periodic orbits
prob = coco_add_slot(prob, 'slot_bd_min_max', @slot_bd_min_max, [], 'bddat');
%prob = coco_add_slot(prob, 'slot_bd_l2', @slot_bd_l2, [], 'bddat');

% Run COCO with the provided continuation problem
bd0 = coco(prob, ...
    'run0', ...             % Name of the continuation run (arbitrary) used for future references
    [], ...                 % Continuation problem is already defined by ode_isol2ep so this term is empty
    'Gamma', ...            % Continuation parameter
    [0.01, 0.25]);          % Parameter range of interest

% Get the data for plotting
p1 = coco_bd_col(bd0, 'Omega'); % Get the Omega parameter
p2 = coco_bd_col(bd0, 'Gamma'); % Get the Gamma parameter
x_max = coco_bd_col(bd0, 'x_max')';   % Get the state vector
stab = coco_bd_col(bd0, 'po.test.USTAB') == 0; % Get the stability

% Plot
plot3(ax, p1, p2, x_max(1, :), 'k-'); % The branch
plot3(ax, p1(stab), p2(stab), x_max(1, stab), 'b.'); % Stable solutions
plot3(ax, p1(~stab), p2(~stab), x_max(1, ~stab), 'r.'); % Unstable solutions
drawnow;

%% Start a continuation for each of the UZ points
uz_labs = coco_bd_labs(bd0, 'UZ');
bd = cell(1, length(uz_labs));
for i = 1:length(uz_labs)

    % Create an empty COCO problem structure
    prob = coco_prob();

    prob = coco_set(prob, 'po', 'TR', false);

    % Create a continuation problem: start with po (periodic orbit) and do
    % a po continuation (periodic orbit), hence use the "po2po" function
    prob = ode_po2po(prob, ...
        '', ...        % OID - only used if embedding this problem in a larger problem
        'run0', ...    % Name of previous run to start from
        uz_labs(i));   % Label on previous run to start from

    % Tell COCO to store extra information about the periodic orbits
    prob = coco_add_slot(prob, 'slot_bd_min_max', @slot_bd_min_max, [], 'bddat');
%     prob = coco_add_slot(prob, 'slot_bd_l2', @slot_bd_l2, [], 'bddat');

    % Turn adaptive meshing on (adapt every n steps)
    prob = coco_set(prob, 'cont', 'NAdapt', 4); % Adapt every 5 steps

    % Increase the maximum step size
    prob = coco_set(prob, 'cont', 'h_max', 10);

    % Increase the number of continuation steps
    prob = coco_set(prob, 'cont', 'ItMX', 300);
    
    
    % Run COCO with the provided continuation problem
    bd{i} = coco(prob, ...
        sprintf('run%d', i), ... % Name of the continuation run (arbitrary) used for future references
        [], ...                  % Continuation problem is already defined by ode_isol2ep so this term is empty
        'Omega', ...             % Continuation parameter
        [0.2, 2]);           % Parameter range of interest

    % Get the data for plotting
    p1 = coco_bd_col(bd{i}, 'Omega'); % Get the Omega parameter
    p2 = coco_bd_col(bd{i}, 'Gamma'); % Get the Gamma parameter
    x_max = coco_bd_col(bd{i}, 'x_max')';   % Get the state vector
    stab = coco_bd_col(bd{i}, 'po.test.USTAB') == 0; % Get the stability

    % Plot
    plot3(ax, p1, p2, x_max(1, :), 'k-'); % The branch
    plot3(ax, p1(stab), p2(stab), x_max(1, stab), 'b.'); % Stable solutions
    plot3(ax, p1(~stab), p2(~stab), x_max(1, ~stab), 'r.'); % Unstable solutions
    drawnow;
    
end

%% Continue a branch of saddle node bifurcations

% Get the label of a saddle node bifurcation point
sn_labs = coco_bd_labs(bd{1}, 'SN');

% Create an empty COCO problem structure
prob = coco_prob();

% Increase the number of mesh points used
prob = coco_set(prob, 'coll', 'NTST', 100);

% Create a continuation problem: start with SN (saddle node bifurcation)
% and do a SN continuation (saddle node bifurcation), hence use the "SN2SN"
% function
prob = ode_SN2SN(prob, ...
    '', ...        % OID - only used if embedding this problem in a larger problem
    'run1', ...    % Name of previous run to start from
    sn_labs(1));   % Label on previous run to start from

% Tell COCO to store extra information about the periodic orbits
prob = coco_add_slot(prob, 'slot_bd_min_max', @slot_bd_min_max, [], 'bddat');
% prob = coco_add_slot(prob, 'slot_bd_l2', @slot_bd_l2, [], 'bddat');

% Turn adaptive meshing on (adapt every n steps)
prob = coco_set(prob, 'cont', 'NAdapt', 4); % Adapt every 5 steps

% Increase the maximum step size
% prob = coco_set(prob, 'cont', 'h_max', 10);

% Increase the number of continuation steps
prob = coco_set(prob, 'cont', 'ItMX', 500);

% Run COCO with the provided continuation problem
bd_sn = coco(prob, ...
    sprintf('run%d', i), ... % Name of the continuation run (arbitrary) used for future references
    [], ...                  % Continuation problem is already defined by ode_isol2ep so this term is empty
    {'Gamma', 'Omega'}, ...  % Continuation parameter
    {[0.01,  0.25], [0.2, 2]}); % Parameter range of interest

% Get the data for plotting
p1 = coco_bd_col(bd_sn, 'Omega'); % Get the Omega parameter
p2 = coco_bd_col(bd_sn, 'Gamma'); % Get the Gamma parameter
x_max = coco_bd_col(bd_sn, 'x_max')';   % Get the state vector

% Plot
plot3(ax, p1, p2, x_max(1, :), 'g.-'); % The branch
drawnow;