% OP_SUM_U_MP: sum up the values of the solution, e.g. for averaging
%
%   sum_u = op_sum_u_mp (spu, msh, patch_list);

% TBD

function sum_u = op_sum_u_mp(spu, msh, u, patch_list)
    
    sum_u = 0;
    
    for iptc = patch_list
        mesh = msh_precompute(msh.msh_patch{iptc});
        sum_u = sum_u + op_Omega(sp_precompute(spu.sp_patch{iptc}, mesh), mesh, u(spu.gnum{iptc}));
    end
end
