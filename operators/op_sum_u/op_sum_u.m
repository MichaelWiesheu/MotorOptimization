% OP_SUM_U: sum up the values of the solution, e.g. for averaging
%
%   sum_u = op_sum_u (spu, msh);

% TBD

function sum_u = op_sum_u(spu, msh, u)
    sum_u = 0;
    jacdet_weights = msh.jacdet .* msh.quad_weights;
   
    for iel = 1:msh.nel
        u_iel = reshape(u(spu.connectivity(:,iel)), 1, 1, []);
        shpu_iel = reshape (spu.shape_functions(:, :, iel), spu.ncomp, msh.nqn, spu.nsh_max);
        jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);

        sum_u = sum_u + sum(sum(shpu_iel.*u_iel, 3).*jacdet_iel, "all");
    end
    
end

