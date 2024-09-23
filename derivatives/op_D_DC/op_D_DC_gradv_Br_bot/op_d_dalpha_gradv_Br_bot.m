
function varargout = op_d_dalpha_gradv_Br_bot(spv, msh, material)

    gradv = spv.shape_function_gradients;

    nu = material.getNuLinear();
    Br = material.Br;
    alpha = material.Angle;
    nbot = [-cos(alpha); -sin(alpha)];

    jacdet = permute(msh.jacdet .* msh.quad_weights, [3,1,4,2]);

    values = reshape(sum(nu.*jacdet.*sum(gradv.*Br.*nbot, 1), 2), [], 1);

    rows = reshape(spv.connectivity, [], 1);
    cols = ones(size(rows));

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, values, spv.ndof, 1);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    else
        error('op_d_dalpha_gradv_Br_bot: wrong number of output arguments')
    end
end