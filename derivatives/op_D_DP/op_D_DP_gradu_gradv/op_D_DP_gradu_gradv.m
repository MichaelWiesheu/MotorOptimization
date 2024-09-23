function varargout = op_D_DP_gradu_gradv(spu, spv, spg, msh, dCdP)

    dKdC = op_D_DC_gradu_gradv(spu, spv, spg, msh);

    dKdP = ttt(dKdC, sptensor(dCdP(:,1:msh.ndim,:)), [3,4], [1,2]);
    
    if (nargout == 1 || nargout == 0)
        varargout{1} = dKdP;
    elseif (nargout == 2)
        [ind, vals] = find(dKdP);
        varargout{1} = ind;
        varargout{2} = vals;
    else
        error ('op_D_DP_gradu_gradv: wrong number of output arguments')
    end
end
