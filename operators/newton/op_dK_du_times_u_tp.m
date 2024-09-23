% OP_DK_DU_TIMES_U_TP: assemble the derivative of the magnetic stiffness matrix A =
% [a(i,j)], a(i,j) = (nu grad u_j, grad v_i) w.r.t magnetic vector
% potential multiplier for the newton scheme exploiting the tensor product structure.
% Directly multiplied by u for faster evaluation
%
%   mat = op_dK_du_times_u_tp(space1, space2, msh, u, material, mu_curve)
%   [indices, values] = op_dK_du_times_u_tp(space1, space2, msh, u, material, mu_curve)
%
% INPUT:
%
%   space1:     object representing the space of trial functions (see sp_vector)
%   space2:     object representing the space of test functions (see sp_vector)
%   msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%   u:       magnetic vector potential multiplier
%   material:   struct containing a B/H curve TBD!
%   mu_curve:   "Iron" so far for TBD!
%
% OUTPUT:
%
%   mat:      assembled stiffness tensor
%   rows:     row indices of the nonzero entries
%   cols:     column indices of the nonzero entries
%   values:   values of the nonzero entries


function varargout = op_dK_du_times_u_tp(space1, space2, msh, u, material)

  for idim = 1:msh.ndim
    size1 = size (space1.sp_univ(idim).connectivity);
    size2 = size (space2.sp_univ(idim).connectivity);
    if ((size1(2) ~= size2(2)) || (size2(2) ~= msh.nel_dir(idim)))
      error ('One of the discrete spaces is not associated to the mesh')
    end
  end
  
  A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col(space1, msh_col, 'value', false, 'gradient', true);
    sp2_col = sp_evaluate_col(space2, msh_col, 'value', false, 'gradient', true);

    A = A + op_dK_du_times_u(sp1_col, sp2_col, msh_col, u, material);
  end

  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 3)
    [rws, cls, vals] = find (A);
    varargout{1} = rws;
    varargout{2} = cls;
    varargout{3} = vals;
  end


end