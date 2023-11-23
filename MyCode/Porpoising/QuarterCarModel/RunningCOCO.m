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

p0 = [180, 50, 30000, 2000, 270000, 500/9, 0.365, 0.365]'; % Initial parameter values
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
    @Suspension, ...     % Right-hand side of equations of motion (optionally the derivatives can also be provided to improve speed)
    t0, ...           % Mesh for the initial solution
    x0, ...           % Initial solution
    {'Ms', 'Mu', 'Ks', 'Cs', 'Kt', 'Vcar', 'A', 'B'}, ... % Parameter names
    p0);              % Initial parameter values

% Add a user defined event to output the values at certain parameter values
prob = coco_add_event(prob, 'UZ', 'Vcar', [100, 150, 200, 250, 300]); % Output when Gamma = 0.05, 0.1 or 0.2

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