%% A helper function to return the L2 norm of a periodic orbit
function [data, res] = slot_bd_l2(prob, data, command, varargin)
% SLOT_BD_MIN_MAX Add the L2 norm to the bifurcation data for each of the
% state variables (individually). It is assumed that the period is one (so
% that a bifurcating branch of periodic orbits has the same L2 norm as the
% underlying Hopf bifurcation point).
%
% Usage: 
%   prob = coco_add_slot(prob, 'slot_bd_l2', @slot_bd_l2, [], 'bddat');
%
% Written by David A.W. Barton (david.barton@bristol.ac.uk) 2015

switch command
    case 'init'
        res = {'x_l2'};
    case 'data'
        chart = varargin{1};
        [data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
        maps = data.coll_seg.maps;
        mesh = data.coll_seg.mesh;
        w = mesh.wts1.*mesh.kas1;
        x = chart.x(uidx(maps.xbp_idx));
        xcn = reshape(maps.W*x, maps.x_shp);
        res  = {sqrt(0.5/maps.NTST*xcn.^2*w')};
    otherwise
        error('Unknown command for slot_bd_l2 - "%s"', command);
end

end