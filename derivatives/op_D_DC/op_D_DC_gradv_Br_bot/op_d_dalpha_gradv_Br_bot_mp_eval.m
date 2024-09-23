
function rhs = op_d_dalpha_gradv_Br_bot_mp_eval(spv, spvEval, mshEval, material, patch_list)

    rhs = zeros (spv.ndof, 1);
    
    for iPatch = patch_list
        rhs_loc = op_d_dalpha_gradv_Br_bot(spvEval{iPatch}, mshEval{iPatch}, material);
    
        if (~isempty (spv.dofs_ornt))
            rhs_loc = spv.dofs_ornt{iptc}(:) .* rhs_loc(:);
        end
    
        rhs(spv.gnum{iPatch}) = rhs(spv.gnum{iPatch}) + rhs_loc;
    
    end
end