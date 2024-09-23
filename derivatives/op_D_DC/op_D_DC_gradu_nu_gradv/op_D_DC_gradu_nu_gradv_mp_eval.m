function A = op_D_DC_gradu_nu_gradv_mp_eval(spu, spuEval, spv, spvEval, spg, spgEval, msh, mshEval, u, material, patch_list)

  inds = cell(msh.npatch, 1);
  vals = cell(msh.npatch, 1);
  
  for iPatch = patch_list
    [ind, v] = op_D_DC_gradu_nu_gradv(spuEval{iPatch}, spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, u(spu.gnum{iPatch}), material);

    ind(:,1) = spv.gnum{iPatch}(ind(:,1));
    ind(:,2) = spu.gnum{iPatch}(ind(:,2));
    ind(:,3) = spg.gnum{iPatch}(ind(:,3));

    if (~isempty (spu.dofs_ornt))
      v = spu.dofs_ornt{iPatch}(ind(:,1))' .* v;
    end
    if (~isempty (spv.dofs_ornt))
      v = spv.dofs_ornt{iPatch}(ind(:,1))' .* v;
    end

    inds{iPatch} = ind;
    vals{iPatch} = v;

  end

  A = sptensor(cell2mat(inds), cell2mat(vals), [spv.ndof, spu.ndof, spg.ndof, msh.ndim]);
end