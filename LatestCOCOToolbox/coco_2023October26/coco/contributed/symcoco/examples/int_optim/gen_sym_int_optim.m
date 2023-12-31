%% Test symbolic differentiation for coco, demo int_optim
% See |coco_folder/po/examples/int_optim| and PO-Tutorial for original |coco| demo, and
% <demo.html> for outputs of demo produced with the derivatives generated
% below.
%
% This demo shows in <demo.html> how one can use symbolic derivatives of up
% to second order for coco computations involving adjoints.
%%
clear
format compact
if sco_isoctave()
    pkg load symbolic
end
%% Create state and parameter names as strings
syms p1 p2 p3 p4         % create symbols for p
syms x1 x2 x3            % create symbols for x
%% ODE
dxdt=[...
    (-p4*(x1^3/3-x1) + (x3-x1)/p2 - x2)/p1;...
    x1-p3;...
    -(x3-x1)/p2];
%% Integrand of objective functional
g=x1/(1+x2^2);
%% Generate code for ODE constraint and integrand of objective functional
sco_sym2funcs(...
    dxdt,...
    {[x1;x2;x3],[p1;p2;p3;p4]},...
    {'x','p'},...
    'filename','sym_mvdP');
sco_sym2funcs(...
    g,...
    {[x1;x2;x3]},...
    {'x'},...
    'filename','sym_g');
