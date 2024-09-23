function varargout = op_D_DC_gradv_Br_bot_tp(space1, space2, msh, material)

    for idim = 1:msh.ndim
        size1 = size (space1.sp_univ(idim).connectivity);
        size2 = size (space2.sp_univ(idim).connectivity);
        if ((size1(2) ~= size2(2)) || (size2(2) ~= msh.nel_dir(idim)))
            error ('One of the discrete spaces is not associated to the mesh')
        end
    end
    
    ind = cell(msh.nel_dir(1), 1);
    vals =  cell(msh.nel_dir(1), 1);
    
    for iel = 1:msh.nel_dir(1)
        msh_col = msh_evaluate_col (msh, iel);
        sp1_col = sp_evaluate_col(space1, msh_col, 'value', true, 'gradient', true);
        sp2_col = sp_evaluate_col(space2, msh_col, 'value', true, 'gradient', true);
    
        [ind{iel}, vals{iel}] = op_D_DC_gradv_Br_bot(sp1_col, sp2_col, msh_col, material);
    end
    
    indices = cell2mat(ind);
    values = cell2mat(vals);

    if (nargout == 1)
        varargout{1} = sptensor(indices, values, [space1.ndof, space2.ndof, msh.ndim]);
    elseif (nargout == 2)
        varargout{1} = indices;
        varargout{2} = values;
    end
end