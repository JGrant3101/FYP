function theme = ep_plot_theme(BT)
%EP_PLOT_THEME    Default plot theme for EP toolbox.
%
%   EP_PLOT_THEME with no arguments displays a list of default themes for
%   all branch types of the EP toolbox.
%
%   EP_PLOT_THEME(BT) returns a struct defining the default plot theme for
%   a branch of type BT. This theme struct contains only settings for EP
%   that override or augment the default COCO plot theme. The full default
%   theme struct is computed in COCO_PLOT_THEME as
%   COCO_MERGE(COCO_PLOT_THEME, EP_PLOT_THEME(BT)).

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: ep_plot_theme.m 3233 2023-05-29 03:05:30Z hdankowicz $

if nargin<1
  display_themes();
else
  
  theme = struct();
  theme.sol.col1 = 'x';
  theme.sol.col2 = 'x';
  theme.sol.RO   = {'bo', 'MarkerFaceColor', 'w', 'MarkerSize', 5};
  theme.plot_sol = @ep_plot_sol;  
  
  switch BT
    
    case 'ep'
      theme.bd.col2  = '||x||_2';
      lspec_s = { 'b-', 'LineWidth', 1, ...
        'DisplayName', 'stable equilibria'};
      lspec_u = {'b--', 'LineWidth', 1, ...
        'DisplayName', 'unstable equilibria'};
      theme.lspec    = {lspec_s, lspec_u};
      theme.ustab    = 'ep.test.USTAB';
      theme.ustabfun = @(x) (x>=1)+1;
      theme.usept    = {'SN', 'HB', 'BP', 'FP'};
      theme.SN       = {'kd', 'MarkerFaceColor', 'g', 'MarkerSize', 8, ...
        'DisplayName', 'Saddle-Node Bifurcation'};
      theme.HB       = {'kd', 'MarkerFaceColor', 'r', 'MarkerSize', 8, ...
        'DisplayName', 'Hopf Bifurcation'};
      theme.NSA      = {'kd', 'MarkerFaceColor', 'm', 'MarkerSize', 8, ...
        'DisplayName', 'Neutral Saddle Point'};     
      theme.sol.SN   = {'go', 'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
        'DisplayName', 'Saddle-Node Bifurcation'};
      theme.sol.HB   = {'ro', 'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
        'DisplayName', 'Hopf Bifurcation'};
      
    case 'ep.SN'
      theme.lspec    = {'g-', 'LineWidth', 1.5, ...
        'DisplayName', 'Saddle-Node Bifurcations'};
      
    case 'ep.HB'
      theme.lspec    = {'r-', 'LineWidth', 1.5, ...
        'DisplayName', 'Hopf Bifurcations'};
      theme.BTP      = {'kp', 'MarkerFaceColor', 'g', 'MarkerSize', 10, ...
        'DisplayName', 'Bogdanov-Takens Point'};
      theme.sol.BTP  = {'go', 'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
        'DisplayName', 'Bogdanov-Takens Point'};
      
    otherwise
      error('%s: unknown solution branch type ''%s''', mfilename, BT);
  end
  
end

end

function display_themes()
tb_info.tb = 'ep';

fprintf('EP default plotting themes:\n');

tb_info.ep.branch_type = 'ep';
default = coco_plot_theme(tb_info);
disp(' ');
disp('''ep'' =');
disp(' ');
disp(default);

tb_info.ep.branch_type = 'ep.SN';
SN = coco_plot_theme(tb_info);
disp(' ');
disp('''ep.SN'' =');
disp(' ');
disp(SN);

tb_info.ep.branch_type = 'ep.HB';
HB = coco_plot_theme(tb_info);
disp(' ');
disp('''ep.HB'' =');
disp(' ');
disp(HB);

end
