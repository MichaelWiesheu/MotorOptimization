function rhs = op_gradu_nu_gradv_times_u_mp_eval(spu, spuEval, spv, spvEval, msh, mshEval, u, material, patch_list)

    rows = cell(msh.npatch, 1);
    cols = cell(msh.npatch, 1);
    vals = cell(msh.npatch, 1);
    
    for iptc = patch_list
        [r, c, v] = op_gradu_nu_gradv_times_u(spuEval{iptc}, spvEval{iptc}, mshEval{iptc}, u(spu.gnum{iptc}), material);
        rows{iptc} = spv.gnum{iptc}(r);
        cols{iptc} = c;
    
        if (~isempty (spv.dofs_ornt))
            v = spv.dofs_ornt{iptc}(r)' .* v;
        end
        if (~isempty (spu.dofs_ornt))
            v = v .* spu.dofs_ornt{iptc}(c)';
        end
        vals{iptc} = v;
    end
    
    rhs = sparse(cell2mat(rows), cell2mat(cols), cell2mat(vals), spv.ndof, 1);
end