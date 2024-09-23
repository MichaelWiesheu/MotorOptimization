% OP_GRADU_NU_GRADV_TP: assemble the stiffness matrix A = [a(i,j)], a(i,j) = (nu grad u_j, grad v_i)
% where nu = nu(B) is given by b-h-curve
%   mat = op_gradu_nu_gradv_tp (spu, spv, msh, u, material, mu_curve);
%   [rows, cols, values] = op_gradu_nu_gradv_tp (spu, spv, msh, u, material, mu_curve);
%
% INPUT:
%
%   spu:     object representing the space of trial functions (see sp_vector)
%   spv:     object representing the space of test functions (see sp_vector)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   u:     magnetic vector potential multiplier
%   material:   struct containing a B/H curve TBD!
%   mu_curve:   "Iron" so far for TBD!
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries

function varargout = op_gradu_nu_gradv_tp(space1, space2, msh, u, material)

    for idim = 1:msh.ndim
        size1 = size (space1.sp_univ(idim).connectivity);
        if (size1(2) ~= msh.nel_dir(idim))
            error ('One of the discrete spaces is not associated to the mesh')
        end
    end
    
    A = spalloc (space1.ndof, space1.ndof, 3*space1.ndof);
    
    for iel = 1:msh.nel_dir(1)
        msh_col = msh_evaluate_col (msh, iel);
        sp1_col = sp_evaluate_col(space1, msh_col, 'value', false, 'gradient', true);
        sp2_col = sp_evaluate_col(space2, msh_col, 'value', false, 'gradient', true);
    
        A = A + op_gradu_nu_gradv(sp1_col, sp2_col, msh_col, u, material);
    end
    
    if (nargout == 1)
        varargout{1} = A;
    elseif (nargout == 3)
        [rows, cols, vals] = find (A);
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = vals;
    end


end