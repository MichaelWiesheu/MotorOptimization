% TBD, faster version for newton solving

function varargout = op_gradu_nu_gradv_times_u(spu, spv, msh, u, material)
    uel = zeros(size(spu.connectivity));
    for iel = 1:msh.nel
        uel(:, iel) = u(spu.connectivity(:, iel));
    end

    gradu = permute(spu.shape_function_gradients, [1,2,5,3,4]); % [dim, nquad, v, u, nel]
    gradv = permute(spv.shape_function_gradients, [1,2,3,5,4]); % [dim, nquad, v, u, nel]

    jacdet = permute(msh.jacdet .* msh.quad_weights, [3,1,4,5,2]);

    gradNj_uj = sum(gradu.*permute(uel, [3,4,5,1,2]), 4);

    B_mag = (sum(gradNj_uj.^2, 1)).^0.5;
    nu = material.getNuNonlinear(B_mag);

    values = reshape(sum(sum(gradNj_uj.*gradv, 1).*jacdet.*nu, 2), [], 1);

    rows = reshape(spv.connectivity, [], 1);
    cols = ones(size(rows));

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, values, spv.ndof, 1);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    else
        error('op_gradu_nu_gradv_times_u: wrong number of output arguments')
    end
end
