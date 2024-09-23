
function varargout = op_D_DP_gradv_Br_bot(spv, spgP, msh, dCdP, material)
    nParams = size(dCdP, 3);
    dCdPel = zeros([size(spgP.connectivity, 1), msh.ndim, size(dCdP, 3), msh.nel]);
    for iel = 1:msh.nel
        dCdPel(:, :, :, iel) = dCdP(spgP.connectivity(:, iel), 1:msh.ndim, :);
    end
    nu = material.getNuLinear();
    alpha = material.Angle;
    Br_bot = material.Br * [-sin(alpha); cos(alpha)];
    
    gradG = permute(spgP.shape_function_gradients, [1,2,5,6,3,4]);
    gradv = permute(spv.shape_function_gradients, [1,2,3,5,6,4]);

    jacdet = permute(msh.jacdet, [3,1,4,5,6,2]);
    weights = permute(msh.quad_weights, [3,1,4,5,6,2]);

    db1dC = -weights.*jacdet.*nu.*sum(gradG.*Br_bot, 1).*gradv;
    db2dC = weights.*jacdet.*nu.*gradG.*sum(Br_bot.*gradv, 1);

    dbdC = sum(db1dC + db2dC, 2);

    dbdP = squeeze(sum(dbdC.* permute(dCdPel, [2,3,5,6,1,4]), [1,5]));

    rows = repmat(permute(spv.connectivity, [3,1,2]), nParams, 1, 1);
    cols = repmat((1:nParams)', 1, size(spv.connectivity, 1), msh.nel);

    dbdP = reshape(dbdP, [], 1);
    rows = reshape(rows, [], 1);
    cols = reshape(cols, [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, dbdP, spv.ndof, nParams);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = dbdP;
    else
        error ('op_D_DP_gradv_Br_bot: wrong number of output arguments')
    end
end
