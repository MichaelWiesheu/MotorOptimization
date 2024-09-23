

function dVdP = op_D_DP_Omega_mp(spg, msh, dCdP, patch_list)
    nParams = size(dCdP, 3);
    dVdP = zeros(nParams, 1);
    
    for iPatch = patch_list   

        values = op_D_DP_Omega_tp(spg.sp_patch{iPatch}, msh.msh_patch{iPatch}, dCdP(spg.gnum{iPatch},:,:));

        dVdP = dVdP + values;
    end
end
