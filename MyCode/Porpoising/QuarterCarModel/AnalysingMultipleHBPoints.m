%% Running Periodic orbit continuation on a variety of HB points.
% In order for this script to work you first need to run the coco_cont
% script
format long
HBlabs = coco_bd_labs(HBbd{1}, 'RO');
HBlabs = unique(HBlabs);
% Take a subsample, every tenth lab
HBlabs = HBlabs(1:10:length(HBlabs));

% Creating empty cell arrays to later put values in
stableparampo = cell(1, length(HBlabs));
unstableparampo = cell(1, length(HBlabs));
stableZsmax = cell(1, length(HBlabs));
unstableZsMax = cell(1, length(HBlabs));
stableZsmin = cell(1, length(HBlabs));
unstableZsmin = cell(1, length(HBlabs));

stableZumax = cell(1, length(HBlabs));
unstableZuMax = cell(1, length(HBlabs));
stableZumin = cell(1, length(HBlabs));
unstableZumin = cell(1, length(HBlabs));

for i = 1:length(HBlabs)
    HBprob = coco_prob();
    HBprob = coco_set(HBprob, 'coll', 'NTST', 200);
    HBprob = coco_set(HBprob, 'po', 'bifus', true);
    HBprob = ode_HB2po(HBprob, '', 'Test1', HBlabs(i));
    HBprob = coco_set(HBprob, 'cont', 'PtMX', 1500);
    HBprob = coco_set(HBprob, 'cont', 'ItMX', 500);
    HBprob = coco_set(HBprob, 'cont', 'h_min', 1e-10);
    HBprob = coco_set(HBprob, 'cont', 'NAdapt', 1, 'h_max', 0.01);

    % Tell COCO to store extra information about the periodic orbits
    HBprob = coco_add_slot(HBprob, 'slot_bd_min_max', @slot_bd_min_max, [], 'bddat');

    HBPo{i} = coco(HBprob, sprintf('HBpo_run%d', i), [], 1, {param 'po.period'}, [0.02 0.07]);

    % Plot PO results for Zs
    bdHB = coco_bd_read(sprintf('HBpo_run%d', i));
    parampo = coco_bd_col(bdHB, param); % Get the current parameter
    x_max = coco_bd_col(bdHB, 'x_max')';   % Get the state vector
    x_min = coco_bd_col(bdHB, 'x_min')';   % Get the state vector
    stabpo = coco_bd_col(bdHB, 'po.test.USTAB') == 0; % Get the stability
    
    % Saving variables
    stableparampo{1, i} = parampo(stabpo)';
    unstableparampo{1, i} = parampo(~stabpo)';
    stableZsmax{1, i} = x_max(1,stabpo)';
    unstableZsMax{1, i} = x_max(1,~stabpo)';
    stableZsmin{1, i} = x_min(1,stabpo)';
    unstableZsmin{1, i} = x_min(1,~stabpo)';

    stableZumax{1, i} = x_max(2,stabpo)';
    unstableZuMax{1, i} = x_max(2,~stabpo)';
    stableZumin{1, i} = x_min(2,stabpo)';
    unstableZumin{1, i} = x_min(2,~stabpo)';

end