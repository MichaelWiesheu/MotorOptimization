% tbd, faster version of op_u_v

function varargout = op_u_v1(spu, spv, msh, func)
    if ~exist("func", "var")
        func = @(x,y,z) ones(size(x));
    end
    if isnumeric(func)
        func = @(x,y,z) func*ones(size(x));
    end
    for idim = 1:msh.rdim
      x{idim} = reshape (msh.geo_map(idim,:,:), msh.nqn, msh.nel);
    end
    coef = permute(func(x{:}), [3,1,4,5,2]);

    shpu = permute(spu.shape_functions, [4,1,5,2,3]); % [dim, nquad, v, u, nel]
    shpv = permute(spv.shape_functions, [4,1,2,5,3]); % [dim, nquad, v, u, nel]

    jacdet = permute(msh.jacdet .* msh.quad_weights, [3,1,4,5,2]);

    values = reshape(sum(sum(shpu.*shpv, 1).*jacdet.*coef, 2), [], 1);

    rows = reshape(repmat(permute(spv.connectivity, [1,3,2]), 1, size(spu.connectivity, 1), 1), [], 1);
    cols = reshape(repmat(permute(spu.connectivity, [3,1,2]), size(spv.connectivity, 1), 1, 1), [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, values, spv.ndof, spu.ndof);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    else
        error('op_gradu_gradv1: wrong number of output arguments')
    end
end
