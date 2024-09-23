% OP_DK_DU_TIMES_U: assemble the derivative of the magnetic stiffness matrix A =
% [a(i, j)], a(i, j) = (nu grad u_j, grad v_i) w.r.t magnetic vector
% potential multiplier for the newton scheme. Directly multiplied by u for
% faster evaluation
%
%   tensor = op_dK_du_times_u (spu, spv, msh, u, material, mu_curve);
%   [indices, values] = op_dK_du_times_u (spu, spv, msh, u, material, mu_curve);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col/sp_precompute)
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   u:     magnetic vector potential multiplier
%   material:   material class with BH curve
%
% OUTPUT:
%
%   mat:    derivative of assembled stiffness matrix magnetic vector potential
%   rows:   row indices the nonzero entries
%   cols:   column indices the nonzero entries
%   values: values of the nonzero entries

function varargout = op_dK_du_times_u(spu, spv, msh, u, material)
    uel = zeros(size(spu.connectivity));
    for iel = 1:msh.nel
        uel(:, iel) = u(spu.connectivity(:, iel));
    end

    gradu = permute(spu.shape_function_gradients, [1,2,5,3,4]); % [dim, nquad, v, u, nel]
    gradv = permute(spv.shape_function_gradients, [1,2,3,5,4]); % [dim, nquad, v, u, nel]

    jacdet = permute(msh.jacdet .* msh.quad_weights, [3,1,4,5,2]);

    gradNi_ui = sum(gradu.*permute(uel, [3,4,5,1,2]), 4);

    B_mag = (sum(gradNi_ui.^2, 1)).^0.5;
    B_cond = 1e-8; % avoid zero division
    B_mag(B_mag<=B_cond) = B_cond;

    dnu_dB = material.getNuPrimeNonlinear(B_mag);
    % nu = material.getNuNonlinear(B_mag);

    values = reshape(sum(sum(gradv.*gradNi_ui, 1).*(sum(gradNi_ui.*gradu, 1)).*jacdet.*dnu_dB./B_mag, 2), [], 1);
    % values = reshape(sum(sum(gradv.*gradu, 1).*nu.*jacdet + sum(gradv.*gradNi_ui, 1).*(sum(gradNi_ui.*gradu, 1)).*jacdet.*dnu_dB./B_mag, 2), [], 1);

    rows = reshape(repmat(permute(spv.connectivity, [1,3,2]), 1, size(spu.connectivity, 1), 1), [], 1);
    cols = reshape(repmat(permute(spu.connectivity, [3,1,2]), size(spv.connectivity, 1), 1, 1), [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, values, spv.ndof, spu.ndof);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    else
        error ('op_dK_du_times_u: wrong number of output arguments')
    end
end

% OLD VERSION WITH ELEMENT LOOP
% function varargout = op_dK_du_times_u(spu, spv, msh, u, material)
% 
%     rows = zeros (msh.nel * spv.nsh_max * spu.nsh_max, 1);
%     cols = zeros (msh.nel * spv.nsh_max * spu.nsh_max, 1);
%     values = zeros (msh.nel * spv.nsh_max * spu.nsh_max, 1);
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
%             B_loc = (sum(gradNi_ui.^2, 1)).^0.5;
%             B_cond = 1e-8; % avoid zero division
%             B_loc(B_loc<=B_cond) = B_cond;
% 
%             dnu_dB = material.getNuPrimeNonlinear(B_loc);
% 
% %             res1 = sum(gradv_iel.*sum(gradu_iel.*u_iel, 4), 1).*dnu_dB.*(Ni_x_ui.*Nj_x + Ni_y_ui.*Nj_y + B_cond)./B_loc;
%             res = sum(gradv_iel.*gradNi_ui, 1).*(sum(gradNi_ui.*gradu_iel, 1));
%             coeff = jacdet_iel.*dnu_dB./B_loc;
% %             res1 = 2*sum(gradv_iel.*sum(gradu_iel.*u_iel1, 4), 1).*dnu_dB2.*(Ni_x_ui.*Nj_x + Ni_y_ui.*Nj_y);
% 
%             elementary_values =  reshape(sum(res.*coeff, 2), spv.nsh_max, spu.nsh_max);
% 
%             [rows_loc, cols_loc] = ndgrid (spv.connectivity(:, iel), spu.connectivity(:, iel));
%             indices = rows_loc & cols_loc;
%             rows(ncounter+(1:spv.nsh(iel)*spu.nsh(iel))) = rows_loc(indices);
%             cols(ncounter+(1:spv.nsh(iel)*spu.nsh(iel))) = cols_loc(indices);
% 
%             values(ncounter+(1:spv.nsh(iel)*spu.nsh(iel))) = elementary_values(indices);
%             ncounter = ncounter + spv.nsh(iel)*spu.nsh(iel);
% 
%         else
%             warning ('geopdes:jacdet_zero_at_quad_node', 'op_dK_du_times_u: singular map in element number %d', iel)
%         end
%     end
% 
%     if (nargout == 1 || nargout == 0)
%         varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
%             values(1:ncounter), spv.ndof, spu.ndof);
%     elseif (nargout == 3)
%         varargout{1} = rows(1:ncounter);
%         varargout{2} = cols(1:ncounter);
%         varargout{3} = values(1:ncounter);
%     else
%         error ('op_dK_du_times_u: wrong number of output arguments')
%     end
% 
% end



