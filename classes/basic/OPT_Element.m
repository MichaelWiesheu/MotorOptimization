classdef OPT_Element < handle
    properties (Access = public)
        OptimizationResponses, NumberOptimizationResponses;
        OptimizationParameters, NumberOptimizationParameters;
        SaveDirectory;
    end

    methods (Abstract)
        f = f_responses(obj, t, x)
        dfdu = partial_f_partial_y(obj, t, x)
        dfdx = partial_f_partial_x(obj, t, x)
        dedx = partial_e_partial_x(obj, t, x)
        optimizationUpdate(obj, x);
    end

    methods
        function obj = OPT_Element()
            obj.resetOptimizationResponses();
            obj.resetOptimizationParameters();
        end

         function addOptimizationParameter(obj, name, initVal, minVal, maxVal)
            pos = numel(obj.OptimizationParameters) + 1;
            if pos > 1
                assert(~any(name == vertcat(obj.OptimizationParameters.Name)), "This parameter has already been added to the element!");
            end
            obj.OptimizationParameters(pos).Name = name;
            obj.OptimizationParameters(pos).Init = initVal;
            obj.OptimizationParameters(pos).Min = minVal;
            obj.OptimizationParameters(pos).Max = maxVal;
            obj.NumberOptimizationParameters = obj.NumberOptimizationParameters + 1;
        end

        function addOptimizationResponse(obj, name)
            pos = numel(obj.OptimizationResponses) + 1;
            if pos > 1
                assert(~any(name == vertcat(obj.OptimizationResponses.Name)), "This response has already been added to the element!");
            end
            obj.OptimizationResponses(pos).Name = name;
            obj.NumberOptimizationResponses = obj.NumberOptimizationResponses + 1;
        end

        function resetOptimizationParameters(obj)
            obj.OptimizationParameters = struct('Name', {}, 'Init', {}, 'Min', {}, 'Max', {});
            obj.NumberOptimizationParameters = 0;
        end

        function resetOptimizationResponses(obj)
            obj.OptimizationResponses = struct('Name', {});
            obj.NumberOptimizationResponses = 0;
        end

    end
end