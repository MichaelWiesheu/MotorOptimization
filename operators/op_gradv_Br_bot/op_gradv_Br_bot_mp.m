
function rhs = op_gradv_Br_bot_mp (space, msh, material, patch_list)

    if (space.npatch ~= msh.npatch)
        error ('op_gradv_Br_bot_mp: the number of patches does not coincide')
    end
    
    rhs = zeros (space.ndof, 1);
    for iptc = patch_list
        rhs_loc = op_gradv_Br_bot_tp (space.sp_patch{iptc}, msh.msh_patch{iptc}, material);
    
        if (~isempty (space.dofs_ornt))
            rhs_loc = space.dofs_ornt{iptc}(:) .* rhs_loc(:);
        end
        rhs(space.gnum{iptc}) = rhs(space.gnum{iptc}) + rhs_loc;
    end
end