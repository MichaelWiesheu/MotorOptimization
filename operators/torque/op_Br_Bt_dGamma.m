% op_Br_Bt_dGamma: calculates the torque resulting from a b_field
%
% INPUT:
%
%   spu:   structure representing the function space (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   coeff: source function evaluated at the quadrature points

function T = op_Br_Bt_dGamma (spu, msh, Az)

    Azel = zeros(size(spu.connectivity));
    for iel = 1:msh.nel
        Azel(:, iel) = Az(spu.connectivity(:, iel));
    end
    Azel = permute(Azel, [3,4,1,2]);

    gradu = spu.shape_function_gradients;

    jacdet = permute(msh.jacdet, [3,1,4,2]);
    weights = permute(msh.quad_weights, [3,1,4,2]);
    normals = permute(msh.normal, [1,2,4,3]);

    r = permute(sum(msh.geo_map.^2, 1).^0.5, [1,2,4,3]);

    % [-By; Bx];
    gradN_Az = sum(gradu.*Azel, 3);
    Bx = gradN_Az(2,:,:,:);
    By = -gradN_Az(1,:,:,:);

    Br = Bx.*normals(1,:,:,:) + By.*normals(2,:,:,:);
    Bt = -Bx.*normals(2,:,:,:) + By.*normals(1,:,:,:);

    T = sum(weights.*jacdet.*r.*Br.*Bt, "all");
end
