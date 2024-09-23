function varargout = op_D_DP_gradu_nu_gradv(spu, spv, spg, msh, dCdP, u, mat)

    dKdC = op_D_DC_gradu_nu_gradv(spu, spv, spg, msh, u, mat);

    dKdP = ttt(dKdC, sptensor(dCdP(:,1:msh.ndim,:)), [3,4], [1,2]);
    
    if (nargout == 1 || nargout == 0)
        varargout{1} = dKdP;
    elseif (nargout == 2)
        [ind, vals] = find(dKdP);
        varargout{1} = ind;
        varargout{2} = vals;
    else
        error ('op_D_DP_gradu_nu_gradv: wrong number of output arguments')
    end
end


% function varargout = op_D_DP_gradu_nu_gradv(spu, spv, spg, msh, u, dCdP, mat)
%     nParams = size(dCdP, 3);
%     uel = zeros(size(spu.connectivity));
%     dCdPel = zeros([size(spg.connectivity, 1), msh.ndim, size(dCdP, 3), msh.nel]);
%     for iel = 1:msh.nel
%         uel(:, iel) = u(spu.connectivity(:, iel));
%         dCdPel(:, :, :, iel) = dCdP(spg.connectivity(:, iel), 1:msh.ndim, :);
%     end
% 
%     gradG = permute(spg.shape_function_gradients, [1,2,5,6,3,4]);
%     gradv = permute(spv.shape_function_gradients, [1,2,3,5,6,4]);
%     gradu = permute(spu.shape_function_gradients, [1,2,5,3,6,4]);
% 
%     jacdet = permute(msh.jacdet, [3,1,4,5,6,2]);
%     weights = permute(msh.quad_weights, [3,1,4,5,6,2]);
% 
%     gradNi_ui = sum(gradu.*permute(uel, [3,4,5,1,6,2]), 4);
%     B_mag = (sum(gradNi_ui.^2, 1)).^0.5;
%     B_cond = 1e-8; % avoid zero division
%     B_mag(B_mag<=B_cond) = B_cond;
%     nu = mat.getNuNonlinear(B_mag);
%     dnudB = mat.getNuPrimeNonlinear(B_mag);
% 
%     dK1dC = -weights.*jacdet.*nu.*sum(gradG.*gradv, 1).*gradu;
%     dK2dC = -weights.*jacdet.*nu.*sum(gradG.*gradu, 1).*gradv;
% 
%     dK3dC = weights.*jacdet.*nu.*gradG.*sum(gradu.*gradv, 1);
% 
%     dK4dC = -weights.*jacdet.*sum(gradu.*gradv, 1).*dnudB./B_mag...
%         .*gradNi_ui.*sum((gradG.*gradNi_ui), 1); 
% 
%     dKdC = sum(dK1dC + dK2dC + dK3dC + dK4dC, 2);
% 
%     dKdP = squeeze(sum(dKdC.* permute(dCdPel, [2,3,5,6,1,4]), [1,5]));
% 
%     rows = repmat(permute(spv.connectivity, [3,1,4,2]), nParams, 1, size(spu.connectivity, 1), 1);
%     cols = repmat(permute(spu.connectivity, [3,4,1,2]), nParams, size(spv.connectivity, 1), 1, 1);
%     dims = repmat((1:nParams)', 1, size(spv.connectivity, 1), size(spu.connectivity, 1), msh.nel);
% 
%     dKdP = reshape(dKdP, [], 1);
%     rows = reshape(rows, [], 1);
%     cols = reshape(cols, [], 1);
%     dims = reshape(dims, [], 1);
% 
%     if (nargout == 1 || nargout == 0)
%         varargout{1} = sptensor ([rows, cols, dims], dKdP, [spv.ndof, spu.ndof, nParams]);
%     elseif (nargout == 2)
%         varargout{1} = [rows, cols, dims];
%         varargout{2} = dKdP;
%     else
%         error ('op_D_DP_gradu_nu_gradv: wrong number of output arguments')
%     end
% end
