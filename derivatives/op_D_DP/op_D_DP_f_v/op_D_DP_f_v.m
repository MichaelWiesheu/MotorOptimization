
function varargout = op_D_DP_f_v(spv, spgP, msh, dCdP)
    nParams = size(dCdP, 3);

    dfdC = full(double(op_D_DC_f_v(spv, spgP, msh)));

    dfdP = squeeze(sum(dfdC.*permute(dCdP(:, 1:msh.ndim, :), [4,1,2,3]), [2,3]));
    [rows, dims, vals] = find(dfdP);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, dims, vals, spv.ndof, nParams);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = dims;
        varargout{3} = vals;
    else
        error ('op_D_DP_f_v: wrong number of output arguments')
    end
end
