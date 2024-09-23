% OP_GRADV_N_BOT_TP: assemble the right-hand side vector r = [r(i)], with  r(i) = (grad(v)_i, n_bot(alpha)).
% where n_bot is the perpendicular normal vector to the magnetization
% direction given by alpha: n_bot = [-sin(alpha); cos(alpha)]
%   rhs = op_d_dalpha_gradv_n_bot_tp (spv, msh, alpha);
% INPUT:

%   spv:   object representing the function space (see sp_vector)
%   msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%   alpha: magnetization angle in rad
%
% OUTPUT:
%
%   rhs:    assembled right hand side

function rhs = op_d_dalpha_gradv_n_bot_tp (space, msh, alpha)

  rhs = zeros (space.ndof, 1);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'gradient', true);

    rhs = rhs + op_d_dalpha_gradv_n_bot (sp_col, msh_col, alpha);
  end

end
