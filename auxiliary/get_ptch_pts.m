% return list of patches and parametric points for physical point

function [sol_patch, pts_phys] = get_ptch_pts(nurbs, x)
    x = cat(1, x, zeros(1, size(x,2)));
    
    pts_phys = zeros(2, size(x,2));
    sol_patch = zeros(1, size(x,2)); %patch in which point is found

    for i = 1:size(x,2)
        for j = 1:numel(nurbs)
            [u_p, converge] = nrbinverse(nurbs(j), x(:,i));
            if (converge)
                sol_patch(i) = j;
                pts_phys(:,i) = u_p;
                break;
            end
        end
    end
    not_found = sum(sol_patch==0);
    if (not_found)
        disp('Could not find ' + string(not_found) + ' points!');
    end
end
