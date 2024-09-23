function varargout = op_D_DP_f_v_tp(space1, space2, msh, dCdP)

    nParams = size(dCdP, 3);

    for idim = 1:msh.ndim
        size1 = size (space1.sp_univ(idim).connectivity);
        size2 = size (space2.sp_univ(idim).connectivity);
        if ((size1(2) ~= size2(2)) || (size2(2) ~= msh.nel_dir(idim)))
            error ('One of the discrete spaces is not associated to the mesh')
        end
    end
    
    rws = cell(msh.nel_dir(1), 1);
    cls = cell(msh.nel_dir(1), 1);
    vals =  cell(msh.nel_dir(1), 1);
    
    for iel = 1:msh.nel_dir(1)
        msh_col = msh_evaluate_col (msh, iel);
        sp1_col = sp_evaluate_col(space1, msh_col, 'value', true, 'gradient', true);
        sp2_col = sp_evaluate_col(space2, msh_col, 'value', true, 'gradient', true);
    
        [rws{iel}, cls{iel}, vals{iel}] = op_D_DP_f_v(sp1_col, sp2_col, msh_col, dCdP);
    end
    
    rows = cell2mat(rws);
    cols = cell2mat(cls);
    values = cell2mat(vals);

    if (nargout == 1)
        varargout{1} = sparse(rows, cols, values, space1.ndof, nParams);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    end
end