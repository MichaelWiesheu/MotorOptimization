classdef IGA_EXCcurrentSine_OPT < IGA_EXCcurrentSine & OPT_Element
    properties (Access = public)
        dRHS
        dCurrent_dIApp;
        dCurrent_dPhase0;
    end

    methods
        function obj = IGA_EXCcurrentSine_OPT(parentElement, phasesStruct, iapp, omega, uvwshift, phase0, nWindings)
           obj@IGA_EXCcurrentSine(parentElement, phasesStruct, iapp, omega, uvwshift, phase0, nWindings);
           obj.dCurrent_dIApp = @(t) sin(obj.Omega*t + obj.UVWshift + obj.Phase0);
           obj.dCurrent_dPhase0 = @(t) obj.IApp*cos(obj.Omega*t + obj.UVWshift + obj.Phase0);
        end
        
        function optimizationUpdate(obj, x)
            if obj.ParentRegion.NumberOptimizationParameters == 0 
                return 
            end
            % Update current
            indCurrent = find(vertcat(obj.ParentRegion.OptimizationParameters.Name) == "IApp");
            if indCurrent > 0
                obj.IApp = x(indCurrent);
            end
            % Update phase
            indPhase0= find(vertcat(obj.ParentRegion.OptimizationParameters.Name) == "Phase0");
            if indPhase0 > 0
                obj.Phase0 = x(indPhase0);
            end
        end

        function dfdy = partial_f_partial_y(obj, t, x, y)
            dfdy = sparse(obj.ParentRegion.NumberOptimizationResponses, obj.ParentRegion.NumberDOF);
        end

        function dedx = partial_e_partial_x(obj, t, x, y)
            % obj.optimizationUpdate(x)
            dedx = sparse(obj.ParentRegion.NumberDOF, obj.ParentRegion.NumberOptimizationParameters);
            
            if obj.ParentRegion.NumberOptimizationParameters == 0 
                return 
            end

            indCurrent = find(vertcat(obj.ParentRegion.OptimizationParameters.Name) == "IApp");
            if indCurrent > 0
                dedx(:, indCurrent) = dedx(:, indCurrent) + obj.dCurrent_dIApp(t)*obj.RHS;
            end
            
            indPhase0 = find(vertcat(obj.ParentRegion.OptimizationParameters.Name) == "Phase0");
            if indPhase0 > 0
                dedx(:, indPhase0) = dedx(:, indPhase0) + obj.dCurrent_dPhase0(t)*obj.RHS;
            end
        end

        function dfdx = partial_f_partial_x(obj, t, x, y)
            dfdx = sparse(obj.ParentRegion.NumberOptimizationResponses, obj.ParentRegion.NumberOptimizationParameters);
        end
        
        function f = f_responses(obj, x)
            f = zeros(obj.ParentRegion.NumberOptimizationResponses, 1);
        end
    end
end