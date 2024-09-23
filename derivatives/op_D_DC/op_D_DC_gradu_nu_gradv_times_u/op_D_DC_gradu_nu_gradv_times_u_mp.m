function dKudC = op_D_DC_gradu_nu_gradv_times_u_mp(spu, spv, spg, msh, u, mat, patch_list)

    indices = cell(msh.npatch, 1);
    values = cell(msh.npatch, 1);
    
    for iPatch = patch_list
        [ind, vals] = op_D_DC_gradu_nu_gradv_times_u_tp(spu.sp_patch{iPatch}, spv.sp_patch{iPatch}, spg.sp_patch{iPatch}, msh.msh_patch{iPatch}, u(spu.gnum{iPatch}), mat);
    
        ind(:,1) = spv.gnum{iPatch}(ind(:,1));
        ind(:,2) = spg.gnum{iPatch}(ind(:,2));
    
        if (~isempty (spv.dofs_ornt))
            vals = spv.dofs_ornt{iPatch}(ind(:,1))' .* vals;
            warning("DOF ornts for u not yet checked");
        end
    
        indices{iPatch} = ind;
        values{iPatch} = vals;
    end
    
    dKudC = sptensor(cell2mat(indices), cell2mat(values), [spv.ndof, spg.ndof, msh.ndim]);
end