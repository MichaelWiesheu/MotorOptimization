function rhs = op_maxwell_force_tp(spu, spv, msh, u, mat)
    
    %   for idim = 1:msh.ndim
    %     size1 = size (spu.sp_univ(idim).connectivity);
    %     size2 = size (spv.sp_univ(idim).connectivity);
    %     if ((size1(2) ~= msh.nel_dir(idim)))
    %       error ('One of the discrete spaces is not associated to the mesh')
    %     end
    %   end
    
    rhs = zeros(spv.ndof, 1);
    
    for iel = 1:msh.nel_dir(1)
        msh_col = msh_evaluate_col (msh, iel);
        spu_col = sp_evaluate_col(spu, msh_col, 'value', true, 'gradient', true);
        spv_col = sp_evaluate_col(spv, msh_col, 'value', true, 'gradient', true);
    
        rhs = rhs + op_maxwell_force(spu_col, spv_col, msh_col, u, mat);
    end
end