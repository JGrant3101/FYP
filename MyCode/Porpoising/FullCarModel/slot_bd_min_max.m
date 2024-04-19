%% A helper function to return the min/max of a periodic orbit
function [data, res] = slot_bd_min_max(prob, data, command, varargin)
% SLOT_BD_MIN_MAX Add the maximum and minimum of the state variables of a
% periodic orbit to the bifurcation data ouput. Note that due to
% discretisation errors the maximum and the minimum are unlikely to vary as
% smoothly as the L2 norm and, as such, may show some anomalous results.
%
% Usage: 
%   prob = coco_add_slot(prob, 'slot_bd_min_max', @slot_bd_min_max, [], 'bddat');
%
% Written by David A.W. Barton (david.barton@bristol.ac.uk) 2015


switch command
    case 'init'
        res = {'x_min', 'x_max'};
    case 'data'
        chart = varargin{1};
        [data, uidx] = coco_get_func_data(prob, 'po.orb.coll', 'data', 'uidx');
        maps = data.coll_seg.maps;
        x = reshape(chart.x(uidx(maps.xbp_idx)), maps.xbp_shp);
        res  = {min(x, [], 2)', max(x, [], 2)'};
    otherwise
        error('Unknown command for slot_bd_min_max - "%s"', command);
end

end
