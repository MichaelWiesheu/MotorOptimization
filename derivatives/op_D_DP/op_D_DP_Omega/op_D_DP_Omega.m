function dVdP = op_D_DP_Omega(spg, msh, dCdP)
    
    dVdC = op_D_DC_Omega(spg, msh);

    dVdP = reshape(sum(full(dVdC).*dCdP(:,1:msh.ndim, :), [1,2]), [], 1);
end

