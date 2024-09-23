function rhs = op_maxwell_force_mp(spu_mp, spv_mp, msh, u, mat, patch_list)
    
    if (spu_mp.npatch ~= msh.npatch)
        error ('op_maxwell_force_mp: the number of patches does not coincide')
    end
    
    rhs = zeros (spv_mp.ndof, 1);
    
    for iptc = patch_list
    
        rhs_loc = op_maxwell_force_tp (spu_mp.sp_patch{iptc}, spv_mp.sp_patch{iptc}, msh.msh_patch{iptc}, u(spu_mp.gnum{iptc}), mat);
    
        if (~isempty (spv_mp.dofs_ornt))
          rhs_loc = spv_mp.dofs_ornt{iptc}(:) .* rhs_loc(:);
        end

        rhs(spv_mp.gnum{iptc}) = rhs(spv_mp.gnum{iptc}) + rhs_loc;
        
    end
end