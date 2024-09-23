function dbdP = op_D_DP_gradv_Br_bot_mp_eval(spv, spvEval, spg, spgEval, msh, mshEval, dCdP, material, patch_list)

    rows = cell(msh.npatch, 1);
    cols = cell(msh.npatch, 1);
    vals = cell(msh.npatch, 1);
    nParams = size(dCdP, 3);
    
    for iPatch = patch_list
        [r, c, v] = op_D_DP_gradv_Br_bot(spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, dCdP(spg.gnum{iPatch},:,:), material);
    
        r = spv.gnum{iPatch}(r);
    
        if (~isempty (spv.dofs_ornt))
            v = spv.dofs_ornt{iPatch}(r)' .* v;
        end
    
        rows{iPatch} = r;
        cols{iPatch} = c;
        vals{iPatch} = v;
    end
    
    dbdP = sparse(cell2mat(rows), cell2mat(cols), cell2mat(vals), spv.ndof, nParams);
end