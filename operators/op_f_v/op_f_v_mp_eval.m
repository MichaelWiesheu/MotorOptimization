
function rhs = op_f_v_mp_eval(spv, spvEval, mshEval, patch_list)

    rhs = zeros (spv.ndof, 1);
    
    for iPatch = patch_list
        % coeff = ones(spvEval{iPatch}.ncomp, mshEval{iPatch}.nqn, mshEval{iPatch}.nel);
        coeff = @(x,y)ones(size(x));

        rhs_loc = op_f_v1(spvEval{iPatch}, mshEval{iPatch}, coeff);
    
        if (~isempty (spv.dofs_ornt))
            rhs_loc = spv.dofs_ornt{iptc}(:) .* rhs_loc(:);
        end
    
        rhs(spv.gnum{iPatch}) = rhs(spv.gnum{iPatch}) + rhs_loc;
    
    end
end