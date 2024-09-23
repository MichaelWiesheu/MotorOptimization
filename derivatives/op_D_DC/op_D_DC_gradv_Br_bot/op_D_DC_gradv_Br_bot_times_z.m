
function varargout = op_D_DC_gradv_Br_bot_times_z(spv, spg, msh, material, z)

    zel = zeros(size(spv.connectivity));
    for iel = 1:msh.nel
        zel(:, iel) = z(spv.connectivity(:, iel));
    end
    zel = permute(zel, [3,4,1,5,6,2]);

    nu = material.getNuLinear();
    alpha = material.Angle;
    Br_bot = material.Br * [-sin(alpha); cos(alpha)];
    
    gradG = permute(spg.shape_function_gradients, [1,2,5,6,3,4]);
    gradNi = permute(spv.shape_function_gradients, [1,2,3,5,6,4]);
    gradNizi = sum(gradNi.*zel, 3);

    jacdet = permute(msh.jacdet, [3,1,4,5,6,2]);
    weights = permute(msh.quad_weights, [3,1,4,5,6,2]);

    db1dC = -weights.*jacdet.*nu.*sum(gradG.*Br_bot, 1).*gradNizi;
    db2dC = weights.*jacdet.*nu.*gradG.*sum(Br_bot.*gradNizi, 1);

    dbdC = squeeze(sum(db1dC + db2dC, 2));

    rows = repmat(permute(spg.connectivity, [3,1,2]), msh.ndim, 1);
    cols = repmat((1:msh.ndim)', 1, spg.nsh_max, msh.nel);

    dbdC = reshape(dbdC, [], 1);
    rows = reshape(rows, [], 1);
    cols = reshape(cols, [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, dbdC, spg.ndof, msh.ndim);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = dbdC;
    else
        error ('op_D_DC_gradv_Br_bot_times_v: wrong number of output arguments')
    end
end
