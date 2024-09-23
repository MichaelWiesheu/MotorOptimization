% MP_EVAL: Compute the value or the derivatives of a function, given by its degrees of freedom, at a given set of points.
%
%   [eu, F] = mp_eval (u, spaces, space_mp, geometry_mp, pts, [options]);
%
function [eu, F_geo] = mp_eval(u, space_mp, geometry_mp, pts, options)
    if (nargin < 5)
        options = {'value'};
    end
    if (~iscell (options))
        options = {options};
    end

    n_mp = numel(geometry_mp);
    eu = cell(n_mp,1);
    F_geo = cell(n_mp,1);

    for i = 1:n_mp
        u_i = u(space_mp.gnum{i});
        pts = {linspace(0,1,geometry_mp(i).nurbs.order(1)*geometry_mp(i).nurbs.number(1)), linspace(0,1,geometry_mp(i).nurbs.order(2)*geometry_mp(i).nurbs.number(2))};
        [eu_i, F_geo_i] = sp_eval(u_i, space_mp.sp_patch{i}, geometry_mp(i), pts, options);
        eu{i} = eu_i; 
        F_geo{i} = F_geo_i;
    end

end



