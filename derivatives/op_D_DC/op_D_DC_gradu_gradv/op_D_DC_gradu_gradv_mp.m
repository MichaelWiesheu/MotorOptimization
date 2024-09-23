
function dKdC = op_D_DC_gradu_gradv_mp(spu, spv, spg, msh, patch_list)
    
    indices = cell(msh.npatch, 1);
    values = cell(msh.npatch, 1);
    
    for iPatch = patch_list
        [ind, vals] = op_D_DC_gradu_gradv_tp(spu.sp_patch{iPatch}, spv.sp_patch{iPatch}, spg.sp_patch{iPatch}, msh.msh_patch{iPatch});
    
        ind(:,1) = spv.gnum{iPatch}(ind(:,1));
        ind(:,2) = spu.gnum{iPatch}(ind(:,2));
        ind(:,3) = spg.gnum{iPatch}(ind(:,3));
    
        if (~isempty (spu.dofs_ornt))
            vals = spu.dofs_ornt{iPatch}(ind(:,1))' .* vals;
        end
        if (~isempty (spv.dofs_ornt))
            vals = spv.dofs_ornt{iPatch}(ind(:,1))' .* vals;
        end
    
        indices{iPatch} = ind;
        values{iPatch} = vals;
    
    end
    
    dKdC = sptensor(cell2mat(indices), cell2mat(values), [spv.ndof, spu.ndof, spg.ndof, msh.ndim]);
end