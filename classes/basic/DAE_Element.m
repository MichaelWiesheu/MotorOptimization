classdef DAE_Element < handle
    properties (Access = public)
        Name;
        NumberDOF;
        % Store the elements, whose internal DOFS are this external DOFS
        ExternalElements = struct('Element', {}, 'LocalIndices', {});
        BoundaryConditions = struct('Boundaries', {}, 'BoundaryCondition', {}, 'dDofs', {}, 'iDofs', {});
        Excitations = struct('Excitation', {}, 'DOFs', {});
        NumberDOFext = 0;
        Solution = struct('Time', {[]}, 'Values', {[]}, 'IsStatic', {[]});
        G, H;           % Boundary Condition matrices
        D, I;           % Dependent and Independent degrees of freedom
    end

    methods (Abstract)
        % Equation system: M(t,y)y' = K(t,y)y + f(t,y)
        MassMatrix(obj, t, y)               % M: Mass matrix
        StiffnessMatrix(obj, t, y)          % K: Stiffness matrix
        ForceVector(obj, t, y)              % f: Right-hand side, which is not in K
        JacobianMatrix(obj, t, y)           % J: Nonlinear contributions of K and RHS
        postprocess(obj, ts, ys)
    end

    methods
        function addExternalElement(obj, element, indices, pos)
            if ~exist("indices", "var")
                indices = 1:element.NumberDOF;
            end
            if ~exist("pos", "var")
                pos = numel(obj.ExternalElements) + 1;
            end
            obj.ExternalElements(pos).Element = element;
            obj.ExternalElements(pos).LocalIndices = indices;
            obj.NumberDOFext = obj.NumberDOFext + numel(indices);
        end

        function resetBoundaryConditions(obj)
            obj.BoundaryConditions = struct('Boundaries', {}, 'BoundaryCondition', {}, 'dDofs', {}, 'iDofs', {});
        end

        function resetExcitations(obj)
            obj.Excitations = struct('Excitation', {}, 'DOFs', {});
        end

        function addBoundaryCondition(obj, boundaryCondition)
            pos = numel(obj.BoundaryConditions)+1;
            boundaryCondition.Name = inputname(2);
            obj.BoundaryConditions(pos).BoundaryCondition = boundaryCondition;
        end

        function addExcitation(obj, excitation)
            pos = numel(obj.Excitations)+1;
            excitation.Name = inputname(2);
            obj.Excitations(pos).Excitation = excitation;
        end

        function [b, db] = getBCvalues(obj, t)
            b = sparse(obj.NumberDOF, 1);
            db = sparse(obj.NumberDOF, 1);
            for iBC = 1:numel(obj.BoundaryConditions)
                BC = obj.BoundaryConditions(iBC).BoundaryCondition;
                localDofsD = BC.getBCindices();
                [b_iel, db_iel] = BC.getBCvalues(t);
                b(localDofsD) = b(localDofsD) + b_iel;
                db(localDofsD) = db(localDofsD) + db_iel;
            end
            dofsD = obj.getBCindices();
            b = b(dofsD);
            db = db(dofsD);
        end

        function [G, H] = getBCmatrices(obj)
            G = sparse(obj.NumberDOF, obj.NumberDOF);       % G = (DxI)
            H = sparse(obj.NumberDOF, obj.NumberDOF);       % H = (DxD)
            for iBC = 1:numel(obj.BoundaryConditions)
                BC = obj.BoundaryConditions(iBC).BoundaryCondition;
                [G_iel, H_iel] = BC.getBCmatrices();
                [localDofsD, localDofsI] = BC.getBCindices();
                G(localDofsD, localDofsI) = G(localDofsD, localDofsI) + G_iel;
                H(localDofsD, localDofsD) = H(localDofsD, localDofsD) + H_iel;
            end
            [dofsD, dofsI] = obj.getBCindices();
            G = G(dofsD, dofsI);
            H = H(dofsD, dofsD);
        end

        function [dofsD, dofsI] = getBCindices(obj)
            if ~isempty(obj.D) || ~isempty(obj.I)
                dofsD = obj.D;
                dofsI = obj.I;
                return
            end
            if isempty(obj.BoundaryConditions)
                dofsD = [];
                dofsI = 1:obj.NumberDOF;
                return
            end
            dofsD = [];
            for iBC = 1:numel(obj.BoundaryConditions)
                BC = obj.BoundaryConditions(iBC).BoundaryCondition;
                localDofsD = BC.getBCindices();
                dofsD = union(dofsD, localDofsD);
            end
            dofsI = setdiff((1:obj.NumberDOF), dofsD);

            obj.D = dofsD;
            obj.I = dofsI;
            assert(numel(dofsD) == numel(unique(dofsD)))
            assert(numel(dofsI) == numel(unique(dofsI)))
            assert(numel(dofsI)+numel(dofsD) == obj.NumberDOF);
        end
    end
end