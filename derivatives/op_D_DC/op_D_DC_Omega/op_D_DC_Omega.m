function varargout = op_D_DC_Omega(spg, msh)

    gradg = spg.shape_function_gradients;

    jacdet = permute(msh.jacdet, [3,1,4,2]);
    weights = permute(msh.quad_weights, [3,1,4,2]);

    dVdC = squeeze(sum(weights.*jacdet.*gradg, 2));

    rows = repmat(permute(spg.connectivity, [3,1,2]), msh.ndim, 1, 1);
    cols = repmat((1:msh.ndim)', 1, size(spg.connectivity, 1), msh.nel);

    dVdC = reshape(dVdC, [], 1);
    rows = reshape(rows, [], 1);
    cols = reshape(cols, [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, dVdC, spg.ndof, msh.ndim);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = dVdC;
    else
        error ('op_D_DC_Omega: wrong number of output arguments')
    end
end

