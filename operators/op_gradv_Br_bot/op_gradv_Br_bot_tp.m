
function rhs = op_gradv_Br_bot_tp (space, msh, material)

    rhs = zeros (space.ndof, 1);
    
    for iel = 1:msh.nel_dir(1)
        msh_col = msh_evaluate_col (msh, iel);
        sp_col  = sp_evaluate_col (space, msh_col, 'gradient', true);
    
        rhs = rhs + op_gradv_Br_bot (sp_col, msh_col, material);
    end

end
