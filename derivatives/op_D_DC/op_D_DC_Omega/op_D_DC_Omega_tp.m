% OP_D_DC_OMEGA_TP: assemble the derivative of the patch volume w.r.t the control points of g.
%
%   mat = op_D_DC_Omega_tp (spg, msh);
%   [row, col, values] = op_D_DC_Omega_tp (spg, msh);
%
% INPUT:
%   spg:     object representing the space of the geometry functions (see sp_scalar)
%   msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   row:    nonzero row entries
%   col:    nonzero column entries
%   values: values of the nonzero entries



function varargout = op_D_DC_Omega_tp(space1, msh)

  for idim = 1:msh.ndim
    size1 = size (space1.sp_univ(idim).connectivity);
    if (size1(2) ~= msh.nel_dir(idim))
      error ('One of the discrete spaces is not associated to the mesh')
    end
  end

  A = sparse(space1.ndof, msh.ndim);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col = sp_evaluate_col(space1, msh_col, 'value', true, 'gradient', true);

    A = A + op_D_DC_Omega(sp_col, msh_col);
  end

  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 3)
    [rows, dims, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = dims;
    varargout{3} = vals;
  end

end