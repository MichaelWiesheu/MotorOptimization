% Faster version of op_f_v with multiple functions for vectorial f output

function varargout = op_fs_v (spv, msh, fs)

    for idim = 1:msh.rdim
        x{idim} = reshape (msh.geo_map(idim, :, :), msh.nqn, msh.nel);
    end    
    fev = fs (x{:});
    nfs = size(fev, 1);

    jacdet = permute(msh.jacdet .* msh.quad_weights, [3,1,4,2]);

    shpv = permute(spv.shape_functions, [4,1,2,3]); % [1, nquad, v, nel]

    values = reshape(sum(jacdet.*shpv.*permute(fev, [1,2,4,3]), 2), [], 1);

    rows = reshape(repmat(permute(spv.connectivity, [3,1,2]), nfs, 1, 1), [], 1);
    cols = reshape(repmat((1:nfs)', 1, size(spv.connectivity, 1), msh.nel), [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, values, spv.ndof, nfs);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    else
        error('op_fs_v: wrong number of output arguments')
    end
end