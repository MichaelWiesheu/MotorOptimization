function dKdP = op_D_DP_gradu_gradv_mp_eval(spu, spuEval, spv, spvEval, spg, spgEval, msh, mshEval, dCdP, patch_list)

  inds = cell(msh.npatch, 1);
  vals = cell(msh.npatch, 1);
  nParams = size(dCdP, 3);
  
  for iPatch = patch_list
    [ind, v] = op_D_DP_gradu_gradv(spuEval{iPatch}, spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, dCdP(spg.gnum{iPatch},:,:));

    ind(:,1) = spv.gnum{iPatch}(ind(:,1));
    ind(:,2) = spu.gnum{iPatch}(ind(:,2));

    if (~isempty (spu.dofs_ornt))
      v = spu.dofs_ornt{iPatch}(ind(:,1))' .* v;
    end
    if (~isempty (spv.dofs_ornt))
      v = spv.dofs_ornt{iPatch}(ind(:,1))' .* v;
    end

    inds{iPatch} = ind;
    vals{iPatch} = v;
  end

  dKdP = sptensor(cell2mat(inds), cell2mat(vals), [spv.ndof, spu.ndof, nParams]);
end