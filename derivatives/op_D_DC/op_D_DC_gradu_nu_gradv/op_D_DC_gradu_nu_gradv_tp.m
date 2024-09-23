%% TBD REWORK with INDICES!

function varargout = op_D_DC_gradu_nu_gradv_tp(space1, space2, space3, msh, u, material)

  for idim = 1:msh.ndim
    size1 = size (space1.sp_univ(idim).connectivity);
    size2 = size (space2.sp_univ(idim).connectivity);
    size3 = size (space3.sp_univ(idim).connectivity);
    if ((size1(2) ~= size2(2)) || (size2(2) ~= size3(2)) || (size3(2) ~= msh.nel_dir(idim)))
      error ('One of the discrete spaces is not associated to the mesh')
    end
  end

  A = sptensor([],[],[space1.ndof, space2.ndof, space3.ndof, msh.ndim]);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col(space1, msh_col, 'value', true, 'gradient', true);
    sp2_col = sp_evaluate_col(space2, msh_col, 'value', true, 'gradient', true);
    sp3_col = sp_evaluate_col(space3, msh_col, 'value', true, 'gradient', true);

    A = A + op_D_DC_gradu_nu_gradv(sp1_col, sp2_col, sp3_col, msh_col, u, material);
  end

  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 2)
    [indices, vals] = find (A);
    varargout{1} = indices;
    varargout{2} = vals;
  end

end