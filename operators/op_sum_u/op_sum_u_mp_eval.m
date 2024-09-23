% OP_SUM_U_MP_EVAL: sum up the values of the solution, e.g. for averaging
%
%   sum_u = op_sum_u_mp_eval (spuEval, mshEval, patch_list);

% TBD

function sum_u = op_sum_u_mp_eval(space, spaceEval, mshEval, u, patch_list)
    
    sum_u = 0;
    
    for iptc = patch_list
        sum_u = sum_u + op_sum_u(spaceEval{iptc}, mshEval{iptc}, u(space.gnum{iptc}));
    end
end
