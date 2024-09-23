% OP_D_DC_GRADU_GRADV: assemble the derivative of the stiffness matrix A = [a(i,j)], a(i,j) = (nu grad u_j, grad v_i) w.r.t the control points of g.
%
%   mat = op_gradu_gradv (spu, spv, spg, msh, u, mat);
%   [indices, values] = op_gradu_gradv (spu, spv, spg, msh, u, mat);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   spg:   structure representing the space of geometry functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   u:     discrete magnetic vector potential
%   mat:   material class with reluctivity functions
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   indices:  indices of the nonzero entries
%   values: values of the nonzero entries

function varargout = op_D_DC_gradu_nu_gradv(spu, spv, spg, msh, u, mat)
    uel = zeros(size(spu.connectivity));
    for iel = 1:msh.nel
        uel(:, iel) = u(spu.connectivity(:, iel));
    end

    gradG = permute(spg.shape_function_gradients, [1,2,5,6,3,4]);
    gradv = permute(spv.shape_function_gradients, [1,2,3,5,6,4]);
    gradu = permute(spu.shape_function_gradients, [1,2,5,3,6,4]);

    % % gradG_iel_mod = reshape (gradG(:,:,:,:,iel), msh.ndim, 1, msh.nqn, 1, 1, spgP.nsh_max);

    jacdet = permute(msh.jacdet, [3,1,4,5,6,2]);
    weights = permute(msh.quad_weights, [3,1,4,5,6,2]);

    gradNi_ui = sum(gradu.*permute(uel, [3,4,5,1,6,2]), 4);
    B_mag = (sum(gradNi_ui.^2, 1)).^0.5;
    B_cond = 1e-8; % avoid zero division
    B_mag(B_mag<=B_cond) = B_cond;
    nu = mat.getNuNonlinear(B_mag);
    dnudB = mat.getNuPrimeNonlinear(B_mag);


    dK1dC = -weights.*jacdet.*nu.*sum(gradG.*gradv, 1).*gradu;
    dK2dC = -weights.*jacdet.*nu.*sum(gradG.*gradu, 1).*gradv;

    dK3dC = weights.*jacdet.*nu.*gradG.*sum(gradu.*gradv, 1);

    dK4dC = -weights.*jacdet.*sum(gradu.*gradv, 1).*dnudB./B_mag...
        .*gradNi_ui.*sum((gradG.*gradNi_ui), 1); 

    values = sum(dK1dC + dK2dC + dK3dC + dK4dC, 2);

    rows = repmat(permute(spv.connectivity, [3,4,1,5,6,2]), msh.ndim, 1, 1, size(spu.connectivity, 1), size(spg.connectivity, 1), 1);
    cols = repmat(permute(spu.connectivity, [3,4,5,1,6,2]), msh.ndim, 1, size(spv.connectivity, 1), 1, size(spg.connectivity, 1), 1);
    tens = repmat(permute(spg.connectivity, [3,4,5,6,1,2]), msh.ndim, 1, size(spv.connectivity, 1), size(spu.connectivity, 1), 1, 1);
    dims = repmat((1:msh.ndim)', 1, 1, size(spv.connectivity, 1), size(spu.connectivity, 1), size(spg.connectivity, 1), msh.nel);


    values = reshape(values, [], 1);
    rows = reshape(rows, [], 1);
    cols = reshape(cols, [], 1);
    tens = reshape(tens, [], 1);
    dims = reshape(dims, [], 1);


    if (nargout == 1 || nargout == 0)
        varargout{1} = sptensor ([rows, cols, tens, dims], values, [spv.ndof, spu.ndof, spg.ndof, msh.ndim]);
    elseif (nargout == 2)
        varargout{1} = [rows, cols, tens, dims];
        varargout{2} = values;
    else
        error ('op_D_DC_gradu_nu_gradv: wrong number of output arguments')
    end
end
