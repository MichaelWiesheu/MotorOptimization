classdef  OPT_System < DAE_System & OPT_Element
    properties
        OptimizationMethod = "interior-point";
        OptimizationObjectives, NumberOptimizationObjectives;
        OptimizationConstraints, NumberOptimizationConstraints;
        LowerBounds, UpperBounds;
        lasty;      % Use last solution as starting values for nonlinear system
        StartTime=tic();
        currentx, currentr; % (Re)use responses for both objective and constraints
        currentxgrad, currentgradr;   
        OptimizationHistory;
    end

    methods
        function obj = OPT_System()
            obj.resetOptimizationObjectives();
            obj.resetOptimizationParameters();
            obj.resetOptimizationConstraints();
        end

        function addOptimizationObjective(obj, name, scaling)
            if ~exist("scaling", "var")
                scaling = 1;
            end
            pos = numel(obj.OptimizationObjectives) + 1;
            obj.OptimizationObjectives(pos).Name = name;
            obj.OptimizationObjectives(pos).Scaling = scaling;
            obj.NumberOptimizationObjectives = obj.NumberOptimizationObjectives + 1;
            % Go through elements, check responses, add Indices to Responses
            obj.updateSystemResponses();
            obj.updateSystemParameters();
        end

        function addOptimizationConstraint(obj, name, type, value, scaling)
            assert(any(type == ["<", ">", "="]));
            if ~exist("scaling", "var")
                scaling = 1;
            end
            pos = numel(obj.OptimizationConstraints) + 1;
            obj.OptimizationConstraints(pos).Name = name;
            obj.OptimizationConstraints(pos).Value = value;
            obj.OptimizationConstraints(pos).Type = type;
            obj.OptimizationConstraints(pos).Scaling = scaling;
            obj.NumberOptimizationConstraints = obj.NumberOptimizationConstraints + 1;
            % Go through elements, check responses, add Indices to Responses
            obj.updateSystemResponses();
            obj.updateSystemParameters();
        end

        function x_opt = optimize(obj)
            obj.OptimizationHistory = struct('iteration', {}, 'funccount', {}, 'fval', {}, 'Responses', {}, 'Objectives', {}, 'Constraints', {}, 'Parameters', {}, 'TimeElapsed', {});
            StartOptTime = string(datetime('now', 'Format', "yyyy-MM-dd-HH-mm"));
            obj.SaveDirectory = "plots/"+StartOptTime;
            mkdir(obj.SaveDirectory);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    obj.Elements(iEl).Element.SaveDirectory = obj.SaveDirectory;
                end
            end
            warning off;

            x0norm= (vertcat(obj.OptimizationParameters.Init) - vertcat(obj.OptimizationParameters.Min))./(vertcat(obj.OptimizationParameters.Max) - vertcat(obj.OptimizationParameters.Min));

            problem.solver = 'fmincon';
            problem.objective = @(x) obj.objective_function(x);
            problem.nonlcon = @(x) obj.constraint_function(x);
            problem.x0 = x0norm;
            problem.options = optimoptions(@fmincon, ... 
                'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
                'Display', 'iter-detailed', ...
                'Algorithm', obj.OptimizationMethod, ...
                'MaxIterations', 300, ...
                'OutputFcn', @(a, b, c) obj.outputFunction(a, b, c), ...
                'PlotFcn', @(a, b, c) obj.plotFunction(a, b, c), ... % 'PlotFcn', @optimplotfval, ...
                'StepTolerance', 1e-8, ... 
                'OptimalityTolerance', 1e-6, ...
                'ConstraintTolerance', 1e-6, ...                
                'ScaleProblem', false, ... % scale manually (otherwise only normalized by initial value)
                'HonorBounds', false);   % do not reset design variables that are at the bounds
            problem.Aeq = [];
            problem.beq = [];
            problem.ub = ones(obj.NumberOptimizationParameters, 1);
            problem.lb = zeros(obj.NumberOptimizationParameters, 1);

            obj.StartTime = tic();
            x_opt = fmincon(problem);

            % Rescale to physical values
            x_opt = x_opt.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;
        end

        % Objective function for optimization
        function [f, gradf] = objective_function(obj, x_norm)
            % Rescale values from 0-1 to real values
            x = x_norm.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;

            if nargout == 2
                if isempty(obj.currentxgrad) || ~all(x == obj.currentxgrad)
                    [r, gradr] = obj.response_function(x);
                    obj.currentx = x;
                    obj.currentxgrad = x;
                    obj.currentr = r;
                    obj.currentgradr = gradr;
                else
                    r = obj.currentr;
                    gradr = obj.currentgradr;
                end
            else    % only calculate response without gradient
                if isempty(obj.currentx) || ~all(x == obj.currentx)
                    r = obj.response_function(x);
                    obj.currentx = x;
                    obj.currentr = r;
                else
                    r = obj.currentr;
                end
            end

            f = 0;
            gradf = zeros(obj.NumberOptimizationParameters, 1);
            for iObjective = 1:numel(obj.OptimizationObjectives)
                ind_obj = obj.OptimizationObjectives(iObjective).Index; 
                scale = obj.OptimizationObjectives(iObjective).Scaling; 
                f = f + r(ind_obj)*scale;

                if nargout == 2
                    gradf = gradf + gradr(ind_obj, :)'*scale;
                end
            end
            % Rescale gradient of f because of parameter normalization
            gradf = gradf .* (obj.UpperBounds-obj.LowerBounds);
        end

        function [c, ceq, gradc, gradceq] = constraint_function(obj, x_norm)
            % Rescale values from 0-1 to real values
            x = x_norm.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;

            if nargout > 2
                if isempty(obj.currentxgrad) || ~all(x == obj.currentxgrad)
                    [r, gradr] = obj.response_function(x);
                    obj.currentx = x;
                    obj.currentxgrad = x;
                    obj.currentr = r;
                    obj.currentgradr = gradr;
                else
                    r = obj.currentr;
                    gradr = obj.currentgradr;
                end
            else    % only calculate response without gradient
                if isempty(obj.currentx) || ~all(x == obj.currentx)
                    r = obj.response_function(x);
                    obj.currentx = x;
                    obj.currentr = r;
                else
                    r = obj.currentr;
                end
            end

            % Make values and derivatives for all constraints, do not
            % distinguish inequality and equality
            allc = zeros(obj.NumberOptimizationConstraints, 1);
            allgradc = zeros(obj.NumberOptimizationParameters, obj.NumberOptimizationConstraints);

            for iConstraint = 1:numel(obj.OptimizationConstraints)
                ind_constr = obj.OptimizationConstraints(iConstraint).Index; 
                scale = obj.OptimizationConstraints(iConstraint).Scaling;
                value = obj.OptimizationConstraints(iConstraint).Value;
                type = obj.OptimizationConstraints(iConstraint).Type;

                if type == "<" || type == "="
                    allc(iConstraint) = (r(ind_constr) - value)*scale;
                    if nargout > 2
                        allgradc(:, iConstraint) = gradr(ind_constr, :)*scale;
                    end
                else
                    allc(iConstraint) = -(r(ind_constr) - value)*scale;
                    if nargout > 2
                        allgradc(:, iConstraint) = -gradr(ind_constr, :)*scale;
                    end
                end
            end
            % Rescale gradient of c because of parameter normalization
            allgradc = allgradc .* (obj.UpperBounds-obj.LowerBounds);
            % get right indices for equality and inequality constraints
            ind_eq = find(vertcat(obj.OptimizationConstraints.Type) == "=");
            ind_ineq = union(find(vertcat(obj.OptimizationConstraints.Type) == ">"), ...
                find(vertcat(obj.OptimizationConstraints.Type) == "<"));
            c = allc(ind_ineq);
            ceq = allc(ind_eq);
            if nargout > 2
                gradc = allgradc(:, ind_ineq);
                gradceq = allgradc(:, ind_eq);
            end
        end
        
        function [r, gradr] = response_function(obj, x, t)
            if ~exist("t", "var")
                t = 0;
            end
            % mytimer("Build", "start");
            obj.optimizationUpdate(x);
            % mytimer("Build", "stop");
            if isempty(obj.currentx) || ~all(x == obj.currentx)
                obj.IsInitialized = false;
            end
            if isempty(obj.lasty)
                y0 = zeros(obj.NumberDOF, 1);
            else
                y0 = obj.lasty;
            end
            % mytimer("Newton", "start");
            [y, Jac] = obj.solveStatic(t, y0, false);
            % mytimer("Newton", "stop");
            obj.lasty = y;
            % mytimer("Responses", "start");
            r = obj.f_responses(t, x, y);
            % mytimer("Responses", "stop");

            if nargout == 2
                % mytimer("Gradient", "start");
                dfdx = obj.partial_f_partial_x(t, x, y);
                dedx = obj.partial_e_partial_x(t, x, y);
                dfdy = obj.partial_f_partial_y(t, x, y);
                % Applying boundary conditions
                dedxBC = dedx(obj.I, :) + obj.G'*(dedx(obj.D, :));
                dfdyBC = dfdy(:, obj.I) + dfdy(:, obj.D)*obj.G;
                % mytimer("Gradient", "stop");
                
                % mytimer("GradientSolve", "start");
                if obj.NumberOptimizationObjectives > obj.NumberOptimizationParameters
                    gradr = dfdx - dfdyBC*(Jac\dedxBC);        % direct method
                else
                    gradr = dfdx - (Jac'\dfdyBC')'*dedxBC;     % adjoint method
                end
                % mytimer("GradientSolve", "stop");
            end
        end

        function updateSystemParameters(obj)
            % Find all defined parameters first
            parameters = [];
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    elParameters = vertcat(obj.Elements(iEl).Element.OptimizationParameters.Name);
                    parameters = unique([parameters; elParameters], 'stable');
                end
            end
            % Find indices for all parameters
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    elParameters = vertcat(obj.Elements(iEl).Element.OptimizationParameters.Name);
                    if isempty(elParameters) 
                        continue 
                    end
                    [~, ~, indices] = intersect(elParameters, parameters, 'stable');
                    obj.Elements(iEl).IndicesParameters = indices';
                end
            end
            % Update System Parameters
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    for iIndex = 1:numel(obj.Elements(iEl).IndicesParameters)
                        index = obj.Elements(iEl).IndicesParameters(iIndex);
                        obj.OptimizationParameters(index).Name = obj.Elements(iEl).Element.OptimizationParameters(iIndex).Name;
                        obj.OptimizationParameters(index).Init = obj.Elements(iEl).Element.OptimizationParameters(iIndex).Init;
                        obj.OptimizationParameters(index).Min = obj.Elements(iEl).Element.OptimizationParameters(iIndex).Min;
                        obj.OptimizationParameters(index).Max = obj.Elements(iEl).Element.OptimizationParameters(iIndex).Max;
                    end
                end
            end
            obj.NumberOptimizationParameters = numel(obj.OptimizationParameters);
            obj.LowerBounds = vertcat(obj.OptimizationParameters.Min);
            obj.UpperBounds = vertcat(obj.OptimizationParameters.Max);
        end

        function updateSystemResponses(obj)
            % Find all defined responses first
            responses = [];
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    elResponses= vertcat(obj.Elements(iEl).Element.OptimizationResponses.Name);
                    responses = unique([responses; elResponses], 'stable');
                end
            end
            % Find indices for all responses
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    elResponses = vertcat(obj.Elements(iEl).Element.OptimizationResponses.Name);
                    if isempty(elResponses)
                        continue
                    end
                    [~, ~, indices] = intersect(elResponses, responses, 'stable');
                    obj.Elements(iEl).IndicesResponses = indices';
                end
            end
            % Update system responses
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    for iIndex = 1:numel(obj.Elements(iEl).IndicesResponses)
                        index = obj.Elements(iEl).IndicesResponses(iIndex);
                        obj.OptimizationResponses(index).Name = obj.Elements(iEl).Element.OptimizationResponses(iIndex).Name;
                    end
                end
            end
            obj.NumberOptimizationResponses = numel(obj.OptimizationResponses);
            % Update system objectives and constraints
            for iObjectives = 1:numel(obj.OptimizationObjectives)
                name = obj.OptimizationObjectives(iObjectives).Name;
                responseNames = vertcat(obj.OptimizationResponses.Name);
                obj.OptimizationObjectives(iObjectives).Index = find(name == responseNames);
            end
            for iConstraint = 1:numel(obj.OptimizationConstraints)
                name = obj.OptimizationConstraints(iConstraint).Name;
                responseNames = vertcat(obj.OptimizationResponses.Name);
                obj.OptimizationConstraints(iConstraint).Index = find(name == responseNames);
            end
        end

        % Update all elements in the optimization 
        function optimizationUpdate(obj, x)
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    obj.Elements(iEl).Element.optimizationUpdate(x(ind_x));
                end
            end
        end

        function f = f_responses(obj, t, x, y)
            f = zeros(obj.NumberOptimizationResponses, 1);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    x_iel = x(ind_x);
                    ind_f = obj.Elements(iEl).IndicesResponses;
                    ind_y = obj.Elements(iEl).Indices;
                    f_iel = obj.Elements(iEl).Element.f_responses(t, x_iel, y_iel);
                    f(ind_f) = f(ind_f) + f_iel;
                end
            end
        end

        function dfdy = partial_f_partial_y(obj, t, x, y)
            dfdy = sparse(obj.NumberOptimizationResponses, obj.NumberDOF);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    ind_y = obj.Elements(iEl).Indices;
                    y_iel = y(ind_y);
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    x_iel = x(ind_x);
                    ind_f = obj.Elements(iEl).IndicesResponses;
                    
                    dfdy_iel = obj.Elements(iEl).Element.partial_f_partial_y(t, x_iel, y_iel);
                    dfdy(ind_f, ind_y) = dfdy(ind_f, ind_y) + dfdy_iel;
                end
            end
        end

        function dfdx = partial_f_partial_x(obj, t, x, y)
            dfdx = sparse(obj.NumberOptimizationResponses, obj.NumberOptimizationParameters);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    ind_f = obj.Elements(iEl).IndicesResponses;
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    x_iel = x(ind_x);
                    dfdx_iel = obj.Elements(iEl).Element.partial_f_partial_x(t, x_iel, y_iel);
                    dfdx(ind_f, ind_x) = dfdx(ind_f, ind_x) + dfdx_iel;
                end
            end
        end

        function dedx = partial_e_partial_x(obj, t, x, y)
            dedx = sparse(obj.NumberDOF, obj.NumberOptimizationParameters);
            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    ind_iel = obj.Elements(iEl).Indices;
                    y_iel = y(ind_iel);
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    x_iel = x(ind_x);
                    dedx_iel = obj.Elements(iEl).Element.partial_e_partial_x(t, x_iel, y_iel);
                    dedx(ind_iel, ind_x) = dedx(ind_iel, ind_x) + dedx_iel;
                end
            end
        end

        function stop = outputFunction(obj, xnorm, optimValues, state)
            stop = false;
            x = xnorm .* (obj.UpperBounds - obj.LowerBounds) + obj.LowerBounds;
            if ~all(x == obj.currentx)
                warning("There might be something wrong here?!")
            end
            r = obj.currentr;
            % Store responses in Output struct
            pos = numel(obj.OptimizationHistory) + 1;
            obj.OptimizationHistory(pos).funccount = optimValues.funccount;
            obj.OptimizationHistory(pos).fval = optimValues.fval;
            obj.OptimizationHistory(pos).iteration = optimValues.iteration;
            obj.OptimizationHistory(pos).Responses = obj.OptimizationResponses;
            for iResponse = 1:obj.NumberOptimizationResponses
                obj.OptimizationHistory(pos).Responses(iResponse).Value = r(iResponse);
            end
            obj.OptimizationHistory(pos).Objectives = obj.OptimizationObjectives;
            for iObjective = 1:obj.NumberOptimizationObjectives
                ind_obj = vertcat(obj.OptimizationResponses.Name) == obj.OptimizationObjectives(iObjective).Name;
                obj.OptimizationHistory(pos).Objectives(iObjective).Value = r(ind_obj);
            end
            obj.OptimizationHistory(pos).Constraints = obj.OptimizationConstraints;
            for iConstraint = 1:obj.NumberOptimizationConstraints
                ind_constr = vertcat(obj.OptimizationResponses.Name) == obj.OptimizationConstraints(iConstraint).Name;
                obj.OptimizationHistory(pos).Constraints(iConstraint).Value = r(ind_constr);
            end
            obj.OptimizationHistory(pos).Parameters = obj.OptimizationParameters;
            for iParameter = 1:obj.NumberOptimizationParameters
                obj.OptimizationHistory(pos).Parameters(iParameter).Value = x(iParameter);
            end
            obj.OptimizationHistory(pos).TimeElapsed = toc(obj.StartTime);
            OptimizationHistory = obj.OptimizationHistory;
            save(obj.SaveDirectory + "/OptimizationHistory.mat", 'OptimizationHistory'); 
        end


        function stop = plotFunction(obj, xnorm, optimValues, state)
            % TBD does scaling make sense here? check if iEl is OPT_System and not scale there?
            x = xnorm.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;

            for iEl = 1:numel(obj.Elements)
                if isa(obj.Elements(iEl).Element, "OPT_Element")
                    ind_x = obj.Elements(iEl).IndicesParameters;
                    x_iel = x(ind_x);
                    obj.Elements(iEl).Element.plotFunction(x_iel, optimValues, state)
                end
            end
            stop = false;
        end

        % Resetting
        function resetOptimizationObjectives(obj)
            obj.OptimizationObjectives = struct('Name', {}, 'Index', {}, 'Scaling', {});
            obj.NumberOptimizationObjectives = 0;
        end
        
        function resetOptimizationConstraints(obj)
            obj.OptimizationConstraints = struct('Name', {}, 'Index', {}, 'Type', {}, 'Value', {}, 'Scaling', {});
            obj.NumberOptimizationConstraints = 0;
        end
    end 
end