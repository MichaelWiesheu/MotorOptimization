% OP_F_V: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i).
%
%   rhs = op_f_v (spv, msh, func);
%
% INPUT:
%
%   spv:   structure representing the function space (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   func:  function handle
%

function varargout = op_f_v1 (spv, msh, func)
    if ~exist("func", "var")
        func = @(x,y,z) ones(size(x));
    end
    for idim = 1:msh.rdim
      x{idim} = reshape (msh.geo_map(idim,:,:), msh.nqn, msh.nel);
    end
    coef = permute(func(x{:}), [1,3,2]);

    shpv  = spv.shape_functions; % [nquad, v, nel]

    jacdet = permute(msh.jacdet .* msh.quad_weights, [1,3,2]);

    values = reshape(sum(jacdet.*coef.*shpv, 1), [], 1);

    rows = reshape(spv.connectivity, [], 1);
    cols = ones(size(rows));

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, values, spv.ndof, 1);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    else
        error('op_f_v1: wrong number of output arguments')
    end
end

