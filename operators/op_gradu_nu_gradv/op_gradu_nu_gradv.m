% OP_GRADU_NU_GRADV: assemble the stiffness matrix A = [a(i, j)], a(i, j) = (nu grad u_j, grad v_i)
% where nu = nu(B) is given by b-h-curve
%   mat = op_gradu_gradv (spu, spv, msh, u, material, mu_curve);
%   [rows, cols, values] = op_gradu_gradv (spu, spv, msh, u, material, mu_curve);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   u:     magnetic vector potential multiplier
%   material:   struct containing a B/H curve TBD!
%   mu_curve:   "Iron" so far for TBD!
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries

function varargout = op_gradu_nu_gradv(spu, spv, msh, u, material)
    uel = zeros(size(spu.connectivity));
    for iel = 1:msh.nel
        uel(:, iel) = u(spu.connectivity(:, iel));
    end

    gradu = permute(spu.shape_function_gradients, [1,2,5,3,4]); % [dim, nquad, v, u, nel]
    gradv = permute(spv.shape_function_gradients, [1,2,3,5,4]); % [dim, nquad, v, u, nel]

    jacdet = permute(msh.jacdet .* msh.quad_weights, [3,1,4,5,2]);

    gradNi_ui = sum(gradu.*permute(uel, [3,4,5,1,2]), 4);

    B_mag = (sum(gradNi_ui.^2, 1)).^0.5;
    nu = material.getNuNonlinear(B_mag);

    values = reshape(sum(sum(gradu.*gradv, 1).*jacdet.*nu, 2), [], 1);

    rows = reshape(repmat(permute(spv.connectivity, [1,3,2]), 1, size(spu.connectivity, 1), 1), [], 1);
    cols = reshape(repmat(permute(spu.connectivity, [3,1,2]), size(spv.connectivity, 1), 1, 1), [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, values, spv.ndof, spu.ndof);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    else
        error('op_gradu_nu_gradv: wrong number of output arguments')
    end
end


%OLD VERSION WITH ELEMENT LOOP
% function varargout = op_gradu_nu_gradv(spu, spv, msh, u, material)
% 
%     rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
%     cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
%     values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
% 
%     jacdet_weights = msh.jacdet .* msh.quad_weights;
% 
%     ncounter = 0;
%     for iel = 1:msh.nel
%         if (all (msh.jacdet(:, iel)))
%             gradu_iel = reshape (spu.shape_function_gradients(:, :, :, iel), msh.ndim, msh.nqn, 1, spu.nsh_max);
%             gradv_iel = reshape (spv.shape_function_gradients(:, :, :, iel), msh.ndim, msh.nqn, spv.nsh_max, 1);
%             u_iel = reshape(u(spu.connectivity(:, iel)), [1, 1, 1, spu.nsh_max]);
%             jacdet_iel = reshape (jacdet_weights(:, iel), [1, msh.nqn, 1, 1]);
% 
%             gradNi_ui = sum(gradu_iel.*u_iel, 4);
% 
%             B_loc = (sum(gradNi_ui.^2, 1)).^0.5;
%             nu = material.getNuNonlinear(B_loc);
% 
%             elementary_values = reshape(sum(sum(gradv_iel.*gradu_iel, 1).*jacdet_iel.*nu, 2), spv.nsh_max, spu.nsh_max);
% 
%             [rows_loc, cols_loc] = ndgrid (spv.connectivity(:, iel), spu.connectivity(:, iel));
%             indices = rows_loc & cols_loc;
%             rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc(indices);
%             cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc(indices);
%             values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = elementary_values(indices);
%             ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);
% 
%         else
%             warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradu_nu_gradv: singular map in element number %d', iel)
%         end
%     end
% 
%     if (nargout == 1 || nargout == 0)
%         varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
%             values(1:ncounter), spu.ndof, spu.ndof);
%     elseif (nargout == 3)
%         varargout{1} = rows(1:ncounter);
%         varargout{2} = cols(1:ncounter);
%         varargout{3} = values(1:ncounter);
%     else
%         error ('op_gradu_nu_gradv: wrong number of output arguments')
%     end
% end