
function varargout = op_D_DC_gradu_gradv_times_u_times_z(spu, spv, spg, msh, u, z)
    uel = zeros(size(spu.connectivity));
    for iel = 1:msh.nel
        uel(:, iel) = u(spu.connectivity(:, iel));
    end

    zel = zeros(size(spv.connectivity));
    for iel = 1:msh.nel
        zel(:, iel) = z(spv.connectivity(:, iel));
    end

    gradG = permute(spg.shape_function_gradients, [1,2,5,6,3,4]);
    gradv = permute(spv.shape_function_gradients, [1,2,3,5,6,4]);
    gradu = permute(spu.shape_function_gradients, [1,2,5,3,6,4]);

    gradu_times_u = sum(gradu.* permute(uel, [3,4,5,1,6,2]), 4);
    gradv_times_v = sum(gradv.* permute(zel, [3,4,1,5,6,2]), 3);

    jacdet = permute(msh.jacdet, [3,1,4,5,6,2]);
    weights = permute(msh.quad_weights, [3,1,4,5,6,2]);

    dK1dC = -weights.*jacdet.*sum(gradG.*gradv_times_v, 1).*gradu_times_u;
    dK2dC = -weights.*jacdet.*sum(gradG.*gradu_times_u, 1).*gradv_times_v;

    dK3dC = weights.*jacdet.*gradG.*sum(gradu_times_u.*gradv_times_v, 1);

    dKdC = squeeze(sum(dK1dC + dK2dC + dK3dC, 2));

    rows = repmat(permute(spg.connectivity, [3,1,2]), msh.ndim, 1, 1);
    cols = repmat((1:msh.ndim)', 1, spg.nsh_max, msh.nel);

    dKdC = reshape(dKdC, [], 1);
    rows = reshape(rows, [], 1);
    cols = reshape(cols, [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse(rows, cols, dKdC, spg.ndof, msh.ndim);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = dKdC;
    else
        error ('op_D_DC_gradu_gradv_times_u_times_v: wrong number of output arguments')
    end
end
