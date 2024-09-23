

function dVdP = op_D_DP_Omega_tp(spaceGeo, msh, dCdP)
    nParams = size(dCdP, 3);
    for idim = 1:msh.ndim
        size1 = size (spaceGeo.sp_univ(idim).connectivity);
        if (size1(2) ~= msh.nel_dir(idim))
            error ('One of the discrete spaces is not associated to the mesh')
        end
    end
    
    dVdP = zeros(nParams, 1);
    
    for iel = 1:msh.nel_dir(1)
        msh_col = msh_evaluate_col (msh, iel);
        sp_col = sp_evaluate_col(spaceGeo, msh_col, 'value', true, 'gradient', true);
    
        dVdP = dVdP + op_D_DP_Omega(sp_col, msh_col, dCdP);
    end
end