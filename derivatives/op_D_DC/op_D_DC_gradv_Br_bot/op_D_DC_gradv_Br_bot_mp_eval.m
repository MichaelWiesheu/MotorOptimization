function dbdC = op_D_DC_gradv_Br_bot_mp_eval(spv, spvEval, spg, spgEval, msh, mshEval, material, patch_list)

    indices = cell(msh.npatch, 1);
    values = cell(msh.npatch, 1);
    
    for iPatch = patch_list
        [ind, vals] = op_D_DC_gradv_Br_bot(spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, material);
    
        ind(:,1) = spv.gnum{iPatch}(ind(:,1));
        ind(:,2) = spg.gnum{iPatch}(ind(:,2));

        if (~isempty (spv.dofs_ornt))
            vals = spv.dofs_ornt{iPatch}(ind(:,1))' .* vals;
        end
    
        indices{iPatch} = ind;
        values{iPatch} = vals;
    end
    
    dbdC = sptensor(cell2mat(indices), cell2mat(values), [spv.ndof, spg.ndof, msh.ndim]);
end