% op_Br_Bt_dGamma: calculates the torque resulting from a b_field
%
% INPUT:
%
%   spv:   structure representing the function space (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   coeff: source function evaluated at the quadrature points


function q = op_gradTn_dGamma (spu, msh, temp)
    tel = zeros(size(spu.connectivity));
    for iel = 1:msh.nel
        tel(:, iel) = temp(spu.connectivity(:, iel));
    end
    tel = permute(tel, [3,4,1,2]);

    gradu = spu.shape_function_gradients;

    graduT = sum(gradu.*tel, 3);

    jacdet = permute(msh.jacdet, [3,1,4,2]);
    weights = permute(msh.quad_weights, [3,1,4,2]);
    normals = permute(msh.normal, [1,2,4,3]);

    q = sum(weights.*jacdet.*graduT.*normals, "all");    
end


