classdef IGA_Region_OPT < OPT_Element & IGA_Region

    properties        
        OptimizationControlPoints;
        GeometryFunction;
        GeometryConstraintFunction; GeometryConstraintName; 
        dCdP;
        currentx;
    end

    methods
        function obj = IGA_Region_OPT()
            
        end

        function optimizationUpdate(obj, x)
            % Check if update is necessary!
            if isempty(obj.currentx)
                % continue with model update
            elseif all(x==obj.currentx)
                % model is already up to date
                return
            end
            % Only rebuild spaces and meshes if 
            if ~isempty(obj.GeometryFunction)

                for iParam = 1:numel(obj.OptimizationParameters)
                    params.(obj.OptimizationParameters(iParam).Name) = x(iParam);
                end
                params.draw_geometry = false;
    
                srf = obj.GeometryFunction(params); % TBD Return also angle of magnets, excitation struct?
    
                % Update surface for optimization control points
                ControlPointsInit = zeros(size(obj.ControlPoints));
                for iPatch  = 1:obj.SpacesGeo.npatch
                    ind_loc = obj.SpacesGeo.gnum{iPatch};
                    ControlPointsInit(ind_loc, 1) = reshape(srf(iPatch).coefs(1, :, :, :)./srf(iPatch).coefs(4, :, :, :), [], 1);
                    ControlPointsInit(ind_loc, 2) = reshape(srf(iPatch).coefs(2, :, :, :)./srf(iPatch).coefs(4, :, :, :), [], 1);
                    ControlPointsInit(ind_loc, 3) = reshape(srf(iPatch).coefs(3, :, :, :)./srf(iPatch).coefs(4, :, :, :), [], 1);
                    ControlPointsInit(ind_loc, 4) = reshape(srf(iPatch).coefs(4, :, :, :), [], 1);
                end
                % Apply control point offsets
                for iCtrlPnt = 1:numel(obj.OptimizationControlPoints)
                    paramNr = find(obj.OptimizationControlPoints(iCtrlPnt).Name == vertcat(obj.OptimizationParameters.Name));
                    dim = obj.OptimizationControlPoints(iCtrlPnt).Dimension;
                    for CtrlPntNr = obj.OptimizationControlPoints(iCtrlPnt).Numbers
                        if dim == "x"
                            ControlPointsInit(CtrlPntNr, 1) = ControlPointsInit(CtrlPntNr, 1) + x(paramNr);
                        elseif dim == "y"
                            ControlPointsInit(CtrlPntNr, 2) = ControlPointsInit(CtrlPntNr, 2) + x(paramNr);
                        elseif dim == "r"
                            angle = atan2(obj.ControlPoints(CtrlPntNr, 2), obj.ControlPoints(CtrlPntNr, 1));
                            ControlPointsInit(CtrlPntNr, 1) = ControlPointsInit(CtrlPntNr, 1) + cos(angle)*x(paramNr);
                            ControlPointsInit(CtrlPntNr, 2) = ControlPointsInit(CtrlPntNr, 2) + sin(angle)*x(paramNr);
                        end
                    end
                end
                % update surface with offsets
                for iPatch  = 1:obj.SpacesGeo.npatch
                    ind_loc = obj.SpacesGeo.gnum{iPatch};
                    srf(iPatch).coefs(1, :, :, :) = reshape(ControlPointsInit(ind_loc, 1).*ControlPointsInit(ind_loc, 4), srf(iPatch).number);
                    srf(iPatch).coefs(2, :, :, :) = reshape(ControlPointsInit(ind_loc, 2).*ControlPointsInit(ind_loc, 4), srf(iPatch).number);
                    srf(iPatch).coefs(3, :, :, :) = reshape(ControlPointsInit(ind_loc, 3).*ControlPointsInit(ind_loc, 4), srf(iPatch).number);
                    srf(iPatch).coefs(4, :, :, :) = reshape(ControlPointsInit(ind_loc, 4), srf(iPatch).number);
                end
    
                obj.updatedCdP(x);
    
                obj.importSurface(srf)  
    
                for iMat = 1:numel(obj.Materials)
                    material = obj.Materials(iMat).Material;
                    for iPatch = 1:numel(obj.Materials(iMat).Patches)
                        patch = obj.Materials(iMat).Patches(iPatch);
                        obj.Geometry(patch).Material = material;
                        obj.Geometry(patch).PlotColor = material.getPlotColor();
                        obj.Geometry(patch).PlotAlpha = material.getPlotAlpha();
                    end
                end
                obj.generateMeshes();
                obj.generateSpaces();
            end

            for iEXC = 1:numel(obj.Excitations)
                if isa(obj.Excitations(iEXC).Excitation, "OPT_Element")
                    obj.Excitations(iEXC).Excitation.optimizationUpdate(x);
                end
            end

            obj.currentx = x;
        end

        function optimizationUpdate2(obj, x)
            obj.optimizationUpdate(x);
            obj.updated2CdP2(x);
        end
        
        function f = f_responses(obj, t, x, y)
            obj.optimizationUpdate(x)

            f = zeros(obj.NumberOptimizationResponses, 1);
            if obj.NumberOptimizationResponses == 0 
                return 
            end
            
            indCost = find(vertcat(obj.OptimizationResponses.Name) == "Cost");
            if indCost > 0
                f(indCost) = obj.Cost();
            end
            indSmoothness= find(vertcat(obj.OptimizationResponses.Name) == "SurfaceSmoothness");
            if indSmoothness > 0
                f(indSmoothness) = obj.surfaceSmoothness(x);
            end

            if ~isempty(obj.GeometryConstraintName)
                xStruct = obj.xToStruct(x);
                c = obj.GeometryConstraintFunction(xStruct);
                indGeometryConstraints = contains(vertcat(obj.OptimizationResponses.Name), obj.GeometryConstraintName);
                f(indGeometryConstraints) = c;
            end
        end
        
        function dfdy = partial_f_partial_y(obj, t, x, y)
            dfdy = sparse(obj.NumberOptimizationResponses, obj.NumberDOF);
        end

        function dfdx = partial_f_partial_x(obj, t, x, y)
            obj.optimizationUpdate(x)

            dfdx = sparse(obj.NumberOptimizationResponses, obj.NumberOptimizationParameters);
            if obj.NumberOptimizationResponses == 0 
                return 
            end
            % Predefined Responses:
            % Cost
            indCost = find(vertcat(obj.OptimizationResponses.Name) == "Cost");
            if indCost > 0
                dfdx(indCost, :) = obj.partial_Cost_partial_x();
            end
            indSmoothness= find(vertcat(obj.OptimizationResponses.Name) == "SurfaceSmoothness");
            if indSmoothness > 0
                [~, dSdx] = obj.surfaceSmoothness(x);
                dfdx(indSmoothness, :) = dSdx;
            end
            % Geometry Constraints: numerical evaluation of gradient
            if ~isempty(obj.GeometryConstraintName)
                indGeometryConstraints = contains(vertcat(obj.OptimizationResponses.Name), obj.GeometryConstraintName);
                dcdx = obj.geometryConstraintGradient(x);
                dfdx(indGeometryConstraints, :) = dcdx;
            end
        end

        function dedx = partial_e_partial_x(obj, t, x, y)
            dedx = sparse(obj.NumberDOF, obj.NumberOptimizationParameters);
        end

        % Specific responses
        function cost = Cost(obj)
            cost = 0;
            for iMat = 1:numel(obj.Materials)
                rho = obj.Materials(iMat).Material.getRho();
                matCost = obj.Materials(iMat).Material.getCost();
                patches = obj.Materials(iMat).Patches;
                cost = cost + obj.Length*rho*matCost*op_Omega_mp_eval(obj.MeshesEval, patches);
            end
        end

        function [p, dSdx] = surfaceSmoothness(obj, x)
            r = zeros(numel(obj.OptimizationControlPoints), 1);
            angles = zeros(numel(obj.OptimizationControlPoints), 1);
            r_ind = zeros(numel(obj.OptimizationControlPoints), 1);
            
            % find all angles and offsets first
            for iCtrlPnt = 1:numel(obj.OptimizationControlPoints)
                CtrlPntNr = obj.OptimizationControlPoints(iCtrlPnt).Numbers(1);
                angles(iCtrlPnt) = atan2(obj.ControlPoints(CtrlPntNr, 2), obj.ControlPoints(CtrlPntNr, 1));
                paramNr = find(obj.OptimizationControlPoints(iCtrlPnt).Name == vertcat(obj.OptimizationParameters.Name));
                r_ind(iCtrlPnt) = paramNr;
                r(iCtrlPnt) = x(paramNr);
            end
            [~, sort_ind] = sort(angles);
            deltar = diff(angles(sort_ind));
            % Penalty matrix which containts angle differences
            K = zeros(numel(r), numel(r));
            for iCtrlPnt = 1:numel(deltar)
                K(iCtrlPnt, iCtrlPnt) = K(iCtrlPnt, iCtrlPnt) + 1./deltar(iCtrlPnt);
                K(iCtrlPnt+1, iCtrlPnt) = K(iCtrlPnt+1, iCtrlPnt) - 1./deltar(iCtrlPnt);
                K(iCtrlPnt, iCtrlPnt+1) = K(iCtrlPnt, iCtrlPnt+1) - 1./deltar(iCtrlPnt);
                K(iCtrlPnt+1, iCtrlPnt+1) = K(iCtrlPnt+1, iCtrlPnt+1) + 1./deltar(iCtrlPnt);
            end
            p = r(sort_ind)'*K*r(sort_ind);

            dSdx = zeros(obj.NumberOptimizationParameters, 1);
            dSdx(r_ind(sort_ind)) = 2*K*r(sort_ind);
        end

        function dCostdP = partial_Cost_partial_x(obj)
            dCostdP = zeros(1, obj.NumberOptimizationParameters);
            for iMat = 1:numel(obj.Materials)
                rho = obj.Materials(iMat).Material.getRho();
                matCost = obj.Materials(iMat).Material.getCost();
                patches = obj.Materials(iMat).Patches;
                dCostdP = dCostdP + obj.Length*rho*matCost*op_D_DP_Omega_mp_eval(obj.SpacesGeo, obj.SpacesGeoEval, obj.MeshesEval, obj.dCdP, patches)';
            end
        end

        function dcdx = geometryConstraintGradient(obj, x)
            xStruct = obj.xToStruct(x);
            c0 = obj.GeometryConstraintFunction(xStruct);
            dcdx = zeros(numel(c0), obj.NumberOptimizationParameters);
            for iParam = 1:obj.NumberOptimizationParameters
                x1 = x;
                Step = (obj.OptimizationParameters(iParam).Max - obj.OptimizationParameters(iParam).Min)*1e-6;
                x1(iParam) = x1(iParam) + Step;
                c1 = obj.GeometryConstraintFunction(obj.xToStruct(x1));
                dcdx(:, iParam) = (c1-c0)/Step;
            end
        end

        function updatedCdP(obj, x)
            % TBD maybe don't make the parameter loop for control points,
            % that do not influente the parameters
            for iParam = 1:numel(obj.OptimizationParameters)
                params.(obj.OptimizationParameters(iParam).Name) = x(iParam);
                % make numerical step size relative to parameter ranges
                steps(iParam) = (obj.OptimizationParameters(iParam).Max - obj.OptimizationParameters(iParam).Min)*1e-6;
            end
            obj.dCdP = getdCdP(obj.GeometryFunction, params, obj.SpacesGeo, steps);
            % check optimization control points
            for iCtrlPnt = 1:numel(obj.OptimizationControlPoints)
                paramNr = find(obj.OptimizationControlPoints(iCtrlPnt).Name == vertcat(obj.OptimizationParameters.Name));
                dim = obj.OptimizationControlPoints(iCtrlPnt).Dimension;
                for CtrlPntNr = obj.OptimizationControlPoints(iCtrlPnt).Numbers
                    if dim == "x"
                        obj.dCdP(CtrlPntNr, 1, paramNr) = 1;
                    elseif dim == "y"
                        obj.dCdP(CtrlPntNr, 2, paramNr) = 1;
                    elseif dim == "r"
                        angle = atan2(obj.ControlPoints(CtrlPntNr, 2), obj.ControlPoints(CtrlPntNr, 1));
                        obj.dCdP(CtrlPntNr, 1, paramNr) = cos(angle);
                        obj.dCdP(CtrlPntNr, 2, paramNr) = sin(angle);
                    end
                end
            end
        end

        function addOptimizationResponse(obj, name)
            addOptimizationResponse@OPT_Element(obj, name);
            obj.currentx = [];
        end

        function addGeometryConstraintFunction(obj, name, func)
            obj.GeometryConstraintFunction = func;
            ndim = numel(func(struct));
            obj.GeometryConstraintName = name;
            for iDim = 1:ndim
                obj.addOptimizationResponse(name+num2str(iDim));
            end
        end

        % Functions for adding responses to the element
        function addOptimizationParameter(obj, name, initVal, minVal, maxVal)
            addOptimizationParameter@OPT_Element(obj, name, initVal, minVal, maxVal);
            obj.currentx = [];
        end

        % dim: "x", "y", "r"
        function addOptimizationControlPointOffset(obj, name, numbers, dim, initVal, minVal, maxVal)
            assert(any(dim==["x", "y", "r"]), "The dimension mus be either x, y or r!");
            obj.addOptimizationParameter(name, initVal, minVal, maxVal);
            pos = numel(obj.OptimizationControlPoints) + 1;
            obj.OptimizationControlPoints(pos).Name = name;
            obj.OptimizationControlPoints(pos).Numbers = numbers;
            obj.OptimizationControlPoints(pos).Dimension = dim;
        end

        function resetOptimizationParameters(obj)
            resetOptimizationParameters@OPT_Element(obj);
            obj.OptimizationControlPoints = struct('Name', {}, 'Numbers', {}, 'Dimension', {});
        end

        function plotOptimizationControlPoints(obj)
            numbers = [];
            for iCtrlPts = 1:numel(obj.OptimizationControlPoints)
                numbers = [numbers, reshape(obj.OptimizationControlPoints(iCtrlPts).Numbers, 1, [])];
            end
            scatter(obj.ControlPoints(numbers, 1), obj.ControlPoints(numbers, 2), "filled", "red")
            text(obj.ControlPoints(numbers, 1), obj.ControlPoints(numbers, 2), string(numbers), "HorizontalAlignment", "center", "VerticalAlignment", "middle");
        end

        function plotFunction(obj, x, optimValues, state)
            obj.optimizationUpdate(x);
            if strcmp(state, 'init')
                obj.FigureGeometry = gcf;
            end
            obj.plotGeometry();
            nrbexport(obj.Surface, obj.SaveDirectory + "/" + obj.Name+"Opt"+string(optimValues.iteration)+".txt")
        end

        function xStruct = xToStruct(obj, x)
            for iParam = 1:obj.NumberOptimizationParameters
                xStruct.(obj.OptimizationParameters(iParam).Name) = x(iParam);
            end
        end
    end
end

 