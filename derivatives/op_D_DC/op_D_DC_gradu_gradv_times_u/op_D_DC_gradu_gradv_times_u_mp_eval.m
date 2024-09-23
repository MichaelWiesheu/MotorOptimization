function dKudC = op_D_DC_gradu_gradv_times_u_mp_eval(spu, spuEval, spv, spvEval, spg, spgEval, msh, mshEval, u, patch_list)

    indices = cell(msh.npatch, 1);
    values = cell(msh.npatch, 1);
    
    for iPatch = patch_list
        [ind, vals] = op_D_DC_gradu_gradv_times_u(spuEval{iPatch}, spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, u(spu.gnum{iPatch}));
    
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