
function varargout = op_D_DC_gradv_Br_bot(spv, spg, msh, material)

    nu = material.getNuLinear();
    alpha = material.Angle;
    Br_bot = material.Br * [-sin(alpha); cos(alpha)];
    
    gradG = permute(spg.shape_function_gradients, [1,2,5,6,3,4]);
    gradv = permute(spv.shape_function_gradients, [1,2,3,5,6,4]);

    jacdet = permute(msh.jacdet, [3,1,4,5,6,2]);
    weights = permute(msh.quad_weights, [3,1,4,5,6,2]);

    db1dC = -weights.*jacdet.*nu.*sum(gradG.*Br_bot, 1).*gradv;
    db2dC = weights.*jacdet.*nu.*gradG.*sum(Br_bot.*gradv, 1);

    dbdC = squeeze(sum(db1dC + db2dC, 2));

    rows = repmat(permute(spv.connectivity, [3,1,4,2]), msh.ndim, 1, spg.nsh_max);
    cols = repmat(permute(spg.connectivity, [3,4,1,2]), msh.ndim, spv.nsh_max, 1);
    dims = repmat((1:msh.ndim)', 1, spv.nsh_max, spg.nsh_max, msh.nel);

    dbdC = reshape(dbdC, [], 1);
    rows = reshape(rows, [], 1);
    cols = reshape(cols, [], 1);
    dims = reshape(dims, [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sptensor ([rows, cols, dims], dbdC, [spv.ndof, spg.ndof, msh.ndim]);
    elseif (nargout == 2)
        varargout{1} = [rows, cols, dims];
        varargout{2} = dbdC;
    else
        error ('op_D_DC_gradv_Br_bot: wrong number of output arguments')
    end
end
