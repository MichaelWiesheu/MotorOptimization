classdef  DAE_System < DAE_Element
    properties
        DAEsolution
        Elements = struct('Element', {}, 'Indices', {}, 'IndicesExt', {}, 'IndicesInt', {});
        Mlin, Klin, Flin;           % Linear contributions of Mass, Stiffness and RHS
        IsInitialized = false;
        tEval, yEval, KEval;
        tSpan;
    end

    methods
        function obj = DAE_System()
            obj.NumberDOF = 0;
        end

        function addElement(obj, element)
            pos = numel(obj.Elements) +1;
            obj.Elements(pos).Element = element;

            obj.Elements(pos).IndicesInt = obj.NumberDOF + (1:element.NumberDOF);
            obj.Elements(pos).IndicesExt = obj.getExternalDOFs(element);
            obj.Elements(pos).Indices = [obj.Elements(pos).IndicesExt, obj.Elements(pos).IndicesInt];

            obj.NumberDOF = obj.NumberDOF + element.NumberDOF;
            obj.IsInitialized = false;
        end

        % TBD calculate BC matrices calculate internal and external DOF indices for all elements
        function initialize(obj)
            if ~obj.IsInitialized
               t0 = 0;
               y0 = zeros(obj.NumberDOF, 1);
               [~, obj.Mlin] = obj.MassMatrix(t0, y0);
               [~, obj.Klin] = obj.StiffnessMatrix(t0, y0);
               [~, obj.Flin] = obj.ForceVector(t0, y0);
               [obj.D, obj.I] = obj.getBCindices();
               [obj.G, obj.H] = obj.getBCmatrices();
            end
            % obj.tEval = t0; 
            % obj.yEval = y0(obj.I); 
            % obj.KEval = K(obj.I, obj.I) + K(obj.I, obj.D)*obj.G + obj.G'*K(obj.D, obj.I) + obj.G'*K(obj.D, obj.D)*obj.G;
            obj.IsInitialized = true;
        end

        function [y0, Jac] = solveStatic(obj, t, y0, verbose)
            obj.initialize();
            if ~exist("t", "var")
                t = 0;
            end
            if ~exist("y0", "var")
                y0 = zeros(obj.NumberDOF, 1);
            end
            if ~exist("verbose", "var")
                verbose = true;
            end
            tol = 1e-6;
            % Init
            % M = obj.MassMatrix(t, y) + obj.Mlin;
            K = obj.StiffnessMatrix(t, y0) + obj.Klin;
            f = obj.ForceVector(t, y0) + obj.Flin;
            % Apply boundary conditions
            [b, db] = obj.getBCvalues(t);
            Kred = K(obj.I, obj.I) + K(obj.I, obj.D)*obj.G + obj.G'*K(obj.D, obj.I) + obj.G'*K(obj.D, obj.D)*obj.G;
            fred = f(obj.I) + obj.G'*f(obj.D) + (K(obj.I, obj.D) + obj.G'*K(obj.D, obj.D))*obj.H*b;% - obj.G'*M(obj.D,obj.D)*obj.H*db;
            y0red = y0(obj.I);
            res = Kred*y0red + fred;
            rnorm = res'*res;
            for iNewton = 1:50
                J = obj.JacobianMatrix(t, y0);
                Jac = Kred + J(obj.I, obj.I) + J(obj.I, obj.D)*obj.G + obj.G'*J(obj.D, obj.I) + obj.G'*J(obj.D, obj.D)*obj.G;
                w = Jac\res;    %search direction

                for iLinesearch = 0:10
                    tau = 0.5^iLinesearch;
                    y1red = y0red - tau*w;
                    % reconstruct
                    y1 = obj.reconstructSolution(t, y1red);
                    % Update system
                    K1 = obj.StiffnessMatrix(t, y1) + obj.Klin;
                    f1 = obj.ForceVector(t, y1) + obj.Flin;
                    % Apply boundary conditions
                    Kred1 = K1(obj.I, obj.I) + K1(obj.I, obj.D)*obj.G + obj.G'*K1(obj.D, obj.I) + obj.G'*K1(obj.D, obj.D)*obj.G;
                    fred1 = f1(obj.I) + obj.G'*f1(obj.D) + (K1(obj.I, obj.D) + obj.G'*K1(obj.D, obj.D))*obj.H*b;% - obj.G'*M(obj.D,obj.D)*obj.H*db;
                    % check residual
                    res1 = Kred1*y1red + fred1;
                    rnorm1 = res1'*res1;
                    if rnorm1<rnorm || rnorm1 < tol
                        Kred = Kred1;
                        res = res1;
                        rnorm = rnorm1;
                        y0 = y1;
                        y0red = y0(obj.I);
                        if verbose
                            disp(['NewtonIter: ' 9 num2str(iNewton) 9 'Residual: ' 9 num2str(rnorm1) 9 'Linesearch ' num2str(iLinesearch)]);
                        end
                        break;
                    end
                end
                if iLinesearch == 10
                    Kred = Kred1;
                    res = res1;
                    rnorm = rnorm1;
                    y0 = y1;
                    y0red = y0(obj.I);
                    if verbose 
                        disp(['NewtonIter: ' 9 num2str(iNewton) 9 'Residual: ' 9 num2str(rnorm1) 9 'Linesearch max value ' num2str(iLinesearch)]);
                    end
                end
                if rnorm1 < tol
                    break
                end
            end

            y0 = obj.reconstructSolution(t, y0red);
            obj.postprocess(t, y0);
        end

        function solveTransient(obj, tSpan, y0, solver, options)
            obj.initialize();
            if ~exist("y0", "var")
                y0 = zeros(obj.NumberDOF, 1);
            elseif isempty(y0)
                y0 = zeros(obj.NumberDOF, 1);
            end
            if ~exist("solver", "var")
                solver = @ode15s;
            elseif isempty(solver)
                solver = @ode15s;
            end
            y0 = y0(obj.I);

            M = obj.massFunction(0, y0);

            standardOptions = odeset('InitialStep', 1e-3, 'RelTol', 1e-3, 'AbsTol', 1e-6, 'Jacobian', @(t, y) obj.jacobianFunction(t, y), ...
                'Mass', M, 'MStateDependence', 'none', 'MassSingular', 'yes', 'Stats', 'on',  'OutputFcn', @(t_, y_, f_) obj.standardOutput(t_, y_, f_));
            if ~exist("solver", "var")
                solver = @ode15s;
            end
            % Overwrite or extend the standard options
            if exist("options", "var")
                field_names = fieldnames (options);
                for iopts  = 1:numel (field_names)
                    if ~ isempty(options.(field_names{iopts}))
                        standardOptions.(field_names{iopts}) = options.(field_names{iopts});
                    end
                end
            end
            obj.tSpan = tSpan;
            obj.DAEsolution = solver(@(t, y) obj.solverFunction(t, y), tSpan, y0, standardOptions);

            ysol = obj.reconstructSolution(obj.DAEsolution.x, obj.DAEsolution.y);
            obj.postprocess(obj.DAEsolution.x, ysol);
        end

        function stop = standardOutput(obj, t, y, flag)
            stop = false;
            switch flag
                case 'init'
                    progressbar('Starting Integration: ')
                case 'done'
                    progressbar('finished')
                otherwise
                    progressbar((t-obj.tSpan(1))/(obj.tSpan(end)-obj.tSpan(1)));
            end
            
        end

        function pst(obj)
            ysol = obj.reconstructSolution(obj.DAEsolution.x, obj.DAEsolution.y);
            obj.postprocess(obj.DAEsolution.x, ysol);
        end

        function res = solverFunction(obj, t, yred)
            y = obj.reconstructSolution(t, yred);
            M = obj.MassMatrix(t, y) + obj.Mlin;
            K = obj.StiffnessMatrix(t, y) + obj.Klin;
            f = obj.ForceVector(t, y) + obj.Flin;
            % Apply boundary conditions
            [b, db] = obj.getBCvalues(t);
            Kred = K(obj.I, obj.I) + K(obj.I, obj.D)*obj.G + obj.G'*K(obj.D, obj.I) + obj.G'*K(obj.D, obj.D)*obj.G;
            fred = f(obj.I) + obj.G'*f(obj.D) + (K(obj.I, obj.D) + obj.G'*K(obj.D, obj.D))*obj.H*b - obj.G'*M(obj.D,obj.D)*obj.H*db;

            obj.tEval= t; obj.yEval=yred; obj.KEval=Kred;

            res = Kred*yred + fred;
        end

        function Mred = massFunction(obj, t, yred)
            y = obj.reconstructSolution(t, yred);
            M = obj.MassMatrix(t, y) + obj.Mlin;
            % Apply boundary conditions
            Mred = M(obj.I, obj.I) + M(obj.I, obj.D)*obj.G + obj.G'*M(obj.D, obj.I) + obj.G'*M(obj.D, obj.D)*obj.G;
        end

        function Kred = stiffnessFunction(obj, t, yred)
            y = obj.reconstructSolution(t, yred);
            % If Stiffness Matrix was recently evaluated
            if isempty(obj.tEval) || isempty(obj.yEval)
                K = obj.StiffnessMatrix(t, y) + obj.Klin;
                Kred = K(obj.I, obj.I) + K(obj.I, obj.D)*obj.G + obj.G'*K(obj.D, obj.I) + obj.G'*K(obj.D, obj.D)*obj.G;
            elseif t == obj.tEval && all(obj.yEval==yred)
                Kred = obj.KEval;
            else
                K = obj.StiffnessMatrix(t, y) + obj.Klin;
                Kred = K(obj.I, obj.I) + K(obj.I, obj.D)*obj.G + obj.G'*K(obj.D, obj.I) + obj.G'*K(obj.D, obj.D)*obj.G;
            end
        end

        function Fred = forceFunction(obj, t, yred)
            y = obj.reconstructSolution(t, yred);
            M = obj.MassMatrix(t, y) + obj.Mlin;
            K = obj.StiffnessMatrix(t, y) + obj.Klin;
            f = obj.ForceVector(t, y) + obj.Flin;
            % Apply boundary conditions
            [b, db] = obj.getBCvalues(t);
            Fred = f(obj.I) + obj.G'*f(obj.D) + (K(obj.I, obj.D) + obj.G'*K(obj.D, obj.D))*obj.H*b - obj.G'*M(obj.D,obj.D)*obj.H*db;
        end

        function Jred = jacobianFunction(obj, t, yred)
            y = obj.reconstructSolution(t, yred);
            Kred = obj.stiffnessFunction(t, yred);
            J = obj.JacobianMatrix(t, y);
            % Apply boundary conditions
            Jred = Kred + J(obj.I, obj.I) + J(obj.I, obj.D)*obj.G + obj.G'*J(obj.D, obj.I) + obj.G'*J(obj.D, obj.D)*obj.G;
        end

        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            if nargout == 1
                Mnonlin = sparse(obj.NumberDOF, obj.NumberDOF);
                for iEl = 1:numel(obj.Elements)
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    Mnonlin_iel = obj.Elements(iEl).Element.MassMatrix(t, y_iel);
                    Mnonlin(ind_iel, ind_iel) = Mnonlin(ind_iel, ind_iel) + Mnonlin_iel;
                end
            else
                Mnonlin = sparse(obj.NumberDOF, obj.NumberDOF);
                Mlin = sparse(obj.NumberDOF, obj.NumberDOF);
                for iEl = 1:numel(obj.Elements)
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    [Mnonlin_iel, Mlin_iel] = obj.Elements(iEl).Element.MassMatrix(t, y_iel);
                    Mnonlin(ind_iel, ind_iel) = Mnonlin(ind_iel, ind_iel) + Mnonlin_iel;
                    Mlin(ind_iel, ind_iel) = Mlin(ind_iel, ind_iel) + Mlin_iel;
                end
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            if nargout == 1
                Knonlin = sparse(obj.NumberDOF, obj.NumberDOF);
                for iEl = 1:numel(obj.Elements)
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    K_iel = obj.Elements(iEl).Element.StiffnessMatrix(t, y_iel);
                    Knonlin(ind_iel, ind_iel) = Knonlin(ind_iel, ind_iel) + K_iel;
                end
            else
                Knonlin = sparse(obj.NumberDOF, obj.NumberDOF);
                Klin = sparse(obj.NumberDOF, obj.NumberDOF);
                for iEl = 1:numel(obj.Elements)
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    [Knonlin_iel, Klin_iel] = obj.Elements(iEl).Element.StiffnessMatrix(t, y_iel);
                    Knonlin(ind_iel, ind_iel) = Knonlin(ind_iel, ind_iel) + Knonlin_iel;
                    Klin(ind_iel, ind_iel) = Klin(ind_iel, ind_iel) + Klin_iel;
                end
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            if nargout == 1
                Fnonlin = sparse(obj.NumberDOF, 1);
                for iEl = 1:numel(obj.Elements)
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    Fnonlin_iel = obj.Elements(iEl).Element.ForceVector(t, y_iel);
                    Fnonlin(ind_iel) = Fnonlin(ind_iel) + Fnonlin_iel;
                end
            else
                Fnonlin = sparse(obj.NumberDOF, 1);
                Flin = sparse(obj.NumberDOF, 1);
                for iEl = 1:numel(obj.Elements)
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    [Fnonlin_iel, Flin_iel] = obj.Elements(iEl).Element.ForceVector(t, y_iel);
                    Fnonlin(ind_iel) = Fnonlin(ind_iel) + Fnonlin_iel;
                    Flin(ind_iel) = Flin(ind_iel) + Flin_iel;
                end
            end
        end

        function J = JacobianMatrix(obj, t, y)
            J = sparse(obj.NumberDOF, obj.NumberDOF);
            for iEl = 1:numel(obj.Elements)
                ind_iel = obj.Elements(iEl).Indices;
                y_iel = y(ind_iel);
                J_iel = obj.Elements(iEl).Element.JacobianMatrix(t, y_iel);
                J(ind_iel, ind_iel) = J(ind_iel, ind_iel) + J_iel;
            end
        end

        function [G, H] = getBCmatrices(obj)
            G = sparse(obj.NumberDOF, obj.NumberDOF);    % G = (DxI)
            H = sparse(obj.NumberDOF, obj.NumberDOF);    % H = (DxD)
            for iEl = 1:numel(obj.Elements)
                checkEl = obj.Elements(iEl).Element;
                [G_iel, H_iel] = checkEl.getBCmatrices();
                [loaclDofsD, localDofsI] = checkEl.getBCindices();
                globalDofsD = obj.Elements(iEl).IndicesInt(loaclDofsD);
                globalDofsI = obj.Elements(iEl).IndicesInt(localDofsI);
                G(globalDofsD, globalDofsI) = G(globalDofsD, globalDofsI) + G_iel;
                H(globalDofsD, globalDofsD) = H(globalDofsD, globalDofsD) + H_iel;
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
            dofsD = [];
            dofsI = [];
            for iEl = 1:numel(obj.Elements)
                checkEl = obj.Elements(iEl).Element;
                globalDofs = obj.Elements(iEl).IndicesInt;
                [localDofsD, localDofsI] = checkEl.getBCindices();
                globalDofsD = globalDofs(localDofsD);
                globalDofsI = globalDofs(localDofsI);
                dofsD = [dofsD, globalDofsD];
                dofsI = [dofsI, globalDofsI];
            end
            obj.D = dofsD;
            obj.I = dofsI;
            assert(numel(dofsD) == numel(unique(dofsD)))
            assert(numel(dofsI) == numel(unique(dofsI)))
        end

        function [b, db] = getBCvalues(obj, t)
            b = sparse(obj.NumberDOF, 1);
            db = sparse(obj.NumberDOF, 1);
            for iEl = 1:numel(obj.Elements)
                checkEl = obj.Elements(iEl).Element;
                [b_iel, db_iel] = checkEl.getBCvalues(t);
                localDofsD = checkEl.getBCindices();
                globalDofsD = obj.Elements(iEl).IndicesInt(localDofsD);
                b(globalDofsD) = b(globalDofsD) + b_iel;
                db(globalDofsD) = db(globalDofsD) + db_iel;
            end
            dofsD = obj.getBCindices();
            b = b(dofsD);
            db = db(dofsD);
        end

        function y = reconstructSolution(obj, ts, yred)
            y = zeros(obj.NumberDOF, numel(ts));
            for ti = 1:numel(ts)
                b = obj.getBCvalues(ts(ti));
                y(obj.I, ti) = yred(:, ti);
                y(obj.D, ti) = obj.G*yred(:, ti) + obj.H*b;
            end
        end

        function extDofs = getExternalDOFs(obj, element)
            extDofs = [];
            for iElExt = 1:numel(element.ExternalElements)
                extEl = element.ExternalElements(iElExt).Element;
                extInd = element.ExternalElements(iElExt).LocalIndices;
                extDofs_iel = obj.getElementDOFs(extEl, extInd);
                if isempty(extDofs_iel)
                    warning("The external element has not yet been added?!! Add the dependencies first!");
                end
                extDofs = [extDofs, extDofs_iel];
            end
            assert(all(extDofs == unique(extDofs, 'stable')));
        end

        function elementIndices = getElementDOFs(obj, element, localIndices)
            elementIndices = [];
            for iEl = 1:numel(obj.Elements)
                checkEl = obj.Elements(iEl).Element;
                globalDofs = obj.Elements(iEl).IndicesInt;
                if isa(checkEl, "DAE_System")
                    if checkEl == element
                        elementIndices = globalDofs(localIndices);
                    else
                        elementIndices = globalDofs(checkEl.getElementDOFs(element, localIndices));
                    end
                elseif isa(checkEl, "DAE_Element")
                    if checkEl == element
                        elementIndices = globalDofs(localIndices);
                    end
                else
                    warning("You need either a DAE_Element or DAE_System here");
                end
                if ~isempty(elementIndices)
                    return
                end
            end
        end

         function idxI = getIndependentIndices(obj, element, localIndices)
             if ~exist("localIndices", "var")
                localIndices = 1:element.NumberDOF;
             end
            elementIndices = obj.getElementDOFs(element, localIndices);
            [~, idxI] = intersect(obj.I, elementIndices);
        end

        function [maxErr, Jac_num, Jac_ana] = checkJacobian(obj, t, y0)
            if ~exist("t", "var")
                t = 0;
            end
            [dofsD, dofsI] = obj.getBCindices();
            if ~exist("y0", "var")
                y0 = zeros(numel(dofsI), 1);
            end
            Step = 1e-8;
            f0 = obj.solverFunction(t, y0);
            Jac_ana = full(obj.jacobianFunction(t, y0));
            Jac_num = zeros(size(Jac_ana));
            progressbar('testing Jac')
            for idof = 1:numel(dofsI)
                y1 = y0;
                y1(idof) = y1(idof) + Step;
                f1 = obj.solverFunction(t, y1);
                Jac_num(:, idof) = (f1-f0)/Step;
                progressbar(idof/numel(dofsI));
            end
            progressbar('done')
            e_rel = 1e-3;
            e_abs = 1e-3;
            err = abs(Jac_num-Jac_ana)./(max(abs(Jac_ana*e_rel), e_abs*ones(size(Jac_ana))));

            [rows, cols] = find(err>1);
            ind = sub2ind(size(Jac_ana), rows, cols);
            a = Jac_ana(ind);
            b = Jac_num(ind);
            maxErr = max(err, [], 'all');
        end

        %%%%%%%% postprocessing

        function postprocess(obj, ts, ys)
            obj.Solution.Time = ts;
            obj.Solution.Values = ys;
            for iEl = 1:numel(obj.Elements)
                obj.Elements(iEl).Element.postprocess(ts, ys(obj.Elements(iEl).Indices, :));
            end
        end
    end
end