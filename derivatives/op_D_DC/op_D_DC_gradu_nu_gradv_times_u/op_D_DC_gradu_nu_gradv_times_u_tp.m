function varargout = op_D_DC_gradu_nu_gradv_times_u_tp(spu, spv, spg, msh, u, mat)

    for idim = 1:msh.ndim
        size1 = size (spu.sp_univ(idim).connectivity);
        size2 = size (spv.sp_univ(idim).connectivity);
        size3 = size (spg.sp_univ(idim).connectivity);
        if ((size1(2) ~= size2(2)) || (size2(2) ~= size3(2)) || (size3(2) ~= msh.nel_dir(idim)))
            error ('One of the discrete spaces is not associated to the mesh')
        end
    end
    
    indices = cell(msh.nel_dir(1), 1);
    values =  cell(msh.nel_dir(1), 1);
    
    for iel = 1:msh.nel_dir(1)
        msh_col = msh_evaluate_col (msh, iel);
        sp1_col = sp_evaluate_col(spu, msh_col, 'value', true, 'gradient', true);
        sp2_col = sp_evaluate_col(spv, msh_col, 'value', true, 'gradient', true);
        sp3_col = sp_evaluate_col(spg, msh_col, 'value', true, 'gradient', true);
    
        [indices{iel}, values{iel}] = op_D_DC_gradu_nu_gradv_times_u(sp1_col, sp2_col, sp3_col, msh_col, u, mat);
    end
    
    indices = cell2mat(indices);
    values = cell2mat(values);

    if (nargout == 1)
        varargout{1} = sptensor(indices, values, [spv.ndof, spg.ndof, msh.ndim]);
    elseif (nargout == 2)
        varargout{1} = indices;
        varargout{2} = values;
    end
end