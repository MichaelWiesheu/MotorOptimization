function A = op_u_v_mp_eval(spu, spuEval, spv, spvEval, msh, mshEval, func, patch_list)

    rows = cell(msh.npatch, 1);
    cols = cell(msh.npatch, 1);
    vals = cell(msh.npatch, 1);
    
    for iPatch = patch_list
        [r, c, v] = op_u_v1(spuEval{iPatch}, spvEval{iPatch}, mshEval{iPatch}, func);
        rows{iPatch} = spv.gnum{iPatch}(r);
        cols{iPatch} = spu.gnum{iPatch}(c);
    
        if (~isempty (spv.dofs_ornt))
            v = spv.dofs_ornt{iPatch}(r)' .* v;
        end
        if (~isempty (spu.dofs_ornt))
            v = v .* spu.dofs_ornt{iPatch}(c)';
        end
        vals{iPatch} = v;
    end
    
    A = sparse(cell2mat(rows), cell2mat(cols), cell2mat(vals), spv.ndof, spu.ndof);

end