
function varargout = op_D_DP_gradu_gradv_times_u(spu, spv, spg, msh, dCdP, u)

    nParams = size(dCdP, 3);
    uel = zeros(size(spu.connectivity));
    dCdPel = zeros([size(spg.connectivity, 1), msh.ndim, size(dCdP, 3), msh.nel]);
    for iel = 1:msh.nel
        uel(:, iel) = u(spu.connectivity(:, iel));
        dCdPel(:, :, :, iel) = dCdP(spg.connectivity(:, iel), 1:msh.ndim, :);
    end

    gradG = permute(spg.shape_function_gradients, [1,2,5,6,3,4]);
    gradv = permute(spv.shape_function_gradients, [1,2,3,5,6,4]);
    gradu = permute(spu.shape_function_gradients, [1,2,5,3,6,4]);

    gradu_times_u = sum(gradu.* permute(uel, [3,4,5,1,6,2]), 4);

    jacdet = permute(msh.jacdet, [3,1,4,5,6,2]);
    weights = permute(msh.quad_weights, [3,1,4,5,6,2]);

    dK1dC = -weights.*jacdet.*sum(gradG.*gradv, 1).*gradu_times_u;
    dK2dC = -weights.*jacdet.*sum(gradG.*gradu_times_u, 1).*gradv;

    dK3dC = weights.*jacdet.*gradG.*sum(gradu_times_u.*gradv, 1);

    dKdC = sum(dK1dC + dK2dC + dK3dC, 2);
    
    dKdP = squeeze(sum(dKdC.* permute(dCdPel, [2,3,5,6,1,4]), [1,5]));

    rows = repmat(permute(spv.connectivity, [3,1,2]), nParams, 1, 1);
    cols = repmat((1:nParams)', 1, size(spv.connectivity, 1), msh.nel);

    dKdP = reshape(dKdP, [], 1);
    rows = reshape(rows, [], 1);
    cols = reshape(cols, [], 1);

    if (nargout == 1 || nargout == 0)
        varargout{1} = sparse (rows, cols, dKdP, spv.ndof, nParams);
    elseif (nargout == 3)
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = dKdP;
    else
        error ('op_D_DP_gradu_gradv_times_u: wrong number of output arguments')
    end

end

%%% More compact but slower version
% function varargout = op_D_DP_gradu_gradv_times_u(spu, spv, spg, msh, u, dCdP)
% 
%     dKdC = op_D_DC_gradu_gradv_times_u(spu, spv, spg, msh, u);
% 
%     dKdP = double(ttt(dKdC, sptensor(dCdP(:,1:msh.ndim,:)), [2,3], [1,2]));
% 
%     if (nargout == 1 || nargout == 0)
%         varargout{1} = dKdP;
%     elseif (nargout == 3)
%         [rows, cols, vals] = find(dKdP);
%         varargout{1} = rows;
%         varargout{2} = cols;
%         varargout{3} = vals;
%     else
%         error ('op_D_DP_gradu_gradv: wrong number of output arguments')
%     end
% end
