

function dVdP = op_D_DP_Omega_mp_eval(spg, spgEval, mshEval, dCdP, patch_list)
    nParams = size(dCdP, 3);
    dVdP = zeros(nParams, 1);
    
    for iPatch = patch_list   

        values = op_D_DP_Omega(spgEval{iPatch}, mshEval{iPatch}, dCdP(spg.gnum{iPatch},:,:));

        dVdP = dVdP + values;
    end
end
