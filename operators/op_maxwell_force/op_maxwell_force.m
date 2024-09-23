
function rhs = op_maxwell_force (spu, spv, msh, u, mat)

    jacdet_weights = msh.jacdet .* msh.quad_weights;
    rhs   = zeros (spv.ndof, 1);
    
    for iel = 1:msh.nel
        if (all (msh.jacdet(:, iel)))
            u_iel = reshape(u(spu.connectivity(:,iel)), [1, 1, spu.nsh_max]);
            jacdet_iel = reshape (jacdet_weights(:,iel), [1, msh.nqn]);
    
            Bx = sum(spu.shape_function_gradients(2,:,:,iel) .* u_iel, 3);
            By = - sum(spu.shape_function_gradients(1,:,:,iel) .* u_iel, 3);
            B = (Bx.^2 + By.^2).^0.5;
    
            if isprop(mat, 'Br')
                Brx = mat.Br*cos(mat.Angle);
                Bry = mat.Br*sin(mat.Angle);
                Hx = (Bx- Brx)/mat.getMuLinear();
                Hy = (By- Bry)/mat.getMuLinear();
                IntBdH = 1/(2*mat.getMuLinear())*(B.^2-Brx^2-Bry^2);
            elseif mat.getType() == "linear"
                Hx = Bx./mat.getMuLinear();
                Hy = By./mat.getMuLinear();
                IntBdH = mat.getCoenergy(B);
            else
                Hx = Bx./mat.getMuNonlinear(B);
                Hy = By./mat.getMuNonlinear(B);
                IntBdH = mat.getCoenergy(B);
            end
            % Maxwell stress tensor
            T(1, 1, :) = (Hx.*Bx - IntBdH);
            T(1, 2, :) = Hx.*By;
            T(2, 1, :) = Bx.*Hy;
            T(2, 2, :) = (Hy.*By - IntBdH);
    
            F = - reshape(sum(T.*spv.shape_function_gradients(:,:,:,:,iel), 2), spv.ncomp, msh.nqn, spv.nsh_max) ;
    
            rhs_loc = sum (sum (F.*jacdet_iel, 1), 2);
    
            indices = find (spv.connectivity(:,iel));
            rhs_loc = rhs_loc(indices); conn_iel = spv.connectivity(indices,iel);
            rhs(conn_iel) = rhs(conn_iel) + rhs_loc(:);
        else
            warning ('geopdes:jacdet_zero_at_quad_node', 'op_d_dp_gradu_gradv: singular map in element number %d', iel)
        end
    end
end



