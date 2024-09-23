
function varargout = op_D_DC_f_v(spv, spg, msh)
    
    shpv = permute(spv.shape_functions, [4,1,2,5,6,3]);

    jacdet = permute(msh.jacdet, [3,1,4,5,6,2]);
    weights = permute(msh.quad_weights, [3,1,4,5,6,2]);

    trDkd = permute(spg.shape_function_gradients, [1,2,5,6,3,4]);

    dfdC = squeeze(sum(weights.*jacdet.*trDkd.*shpv, 2));

    rows = repmat(permute(spv.connectivity, [3,1,4,2]), msh.ndim, 1, spg.nsh_max);
    cols = repmat(permute(spg.connectivity, [3,4,1,2]), msh.ndim, spv.nsh_max, 1);
    dims = repmat((1:msh.ndim)', 1, spv.nsh_max, spg.nsh_max, msh.nel);

    dfdC = reshape(dfdC, [], 1);
    rows = reshape(rows, [], 1);
    cols = reshape(cols, [], 1);
    dims = reshape(dims, [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sptensor ([rows, cols, dims], dfdC, [spv.ndof, spg.ndof, msh.ndim]);
    elseif (nargout == 2)
        varargout{1} = [rows, cols, dims];
        varargout{2} = dfdC;
    else
        error ('op_D_DC_f_v: wrong number of output arguments')
    end
end
