function varargout = op_D_DC_gradu_gradv(spu, spv, spg, msh)

    gradG = permute(spg.shape_function_gradients, [1,2,5,6,3,4]);
    gradv = permute(spv.shape_function_gradients, [1,2,3,5,6,4]);
    gradu = permute(spu.shape_function_gradients, [1,2,5,3,6,4]);

    jacdet = permute(msh.jacdet, [3,1,4,5,6,2]);
    weights = permute(msh.quad_weights, [3,1,4,5,6,2]);

    dK1dC = -weights.*jacdet.*sum(gradG.*gradv, 1).*gradu;
    dK2dC = -weights.*jacdet.*sum(gradG.*gradu, 1).*gradv;

    dK3dC = weights.*jacdet.*gradG.*sum(gradu.*gradv, 1);

    dKdC = squeeze(sum(dK1dC + dK2dC + dK3dC, 2));


    rows = repmat(permute(spv.connectivity, [3,1,4,5,2]), msh.ndim, 1, spu.nsh_max, spg.nsh_max, 1);
    cols = repmat(permute(spu.connectivity, [3,4,1,5,2]), msh.ndim, spv.nsh_max, 1, spg.nsh_max, 1);
    tens = repmat(permute(spg.connectivity, [3,4,5,1,2]), msh.ndim, spv.nsh_max, spu.nsh_max, 1, 1);
    dims = repmat((1:msh.ndim)', 1, spv.nsh_max, spu.nsh_max, spg.nsh_max, msh.nel);

    dKdC = reshape(dKdC, [], 1);
    rows = reshape(rows, [], 1);
    cols = reshape(cols, [], 1);
    tens = reshape(tens, [], 1);
    dims = reshape(dims, [], 1);


    if (nargout == 1 || nargout == 0)
        varargout{1} = sptensor ([rows, cols, tens, dims], dKdC, [spv.ndof, spu.ndof, spg.ndof, msh.ndim]);
    elseif (nargout == 2)
        varargout{1} = [rows, cols, tens, dims];
        varargout{2} = dKdC;
    else
        error ('op_D_DP_gradu_gradv: wrong number of output arguments')
    end
end
