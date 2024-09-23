classdef IGA_MotorClass_OPT < IGA_MotorClass & OPT_System

    properties
        OptimizationAngles;
        LastTorques;
        LastTorquesx;
    end

    methods
        function obj = IGA_MotorClass_OPT(rotor, stator)
            obj@IGA_MotorClass(rotor, stator);
        end

        function [r, gradr] = response_function(obj, x)

            indTorque = find(vertcat(obj.OptimizationResponses.Name) == "Torque");
            assert(indTorque>0, "Torque needs to be a response in order to calculate mean and std!");
            indTorqueMean = find(vertcat(obj.OptimizationResponses.Name) == "TorqueMEAN");
            indTorqueStd = find(vertcat(obj.OptimizationResponses.Name) == "TorqueSTD");
            
            Torques = zeros(numel(obj.OptimizationAngles), 1);
            dTorquesdx = zeros(numel(obj.OptimizationAngles), obj.NumberOptimizationParameters);

            for iAngle = 1:numel(obj.OptimizationAngles)
                theta = obj.OptimizationAngles(iAngle);
                obj.Coupling.setRotationAngle(theta)
                if nargout == 2
                    [ri, gradri] = response_function@OPT_System(obj, x, theta);
                    Torques(iAngle) = ri(indTorque);
                    dTorquesdx(iAngle, :) = gradri(indTorque, :);
                    % Use values of first rotation angle for all other
                    % responses which are not mean or std torque
                    if iAngle == 1
                        r = ri;
                        gradr = gradri;
                    end
                else
                    ri = response_function@OPT_System(obj, x, theta);
                    Torques(iAngle) = ri(indTorque);
                    if iAngle == 1
                        r = ri;
                    end
                end
            end
            Tmean = mean(Torques);
            Tstd = std(Torques, 1, 1);
            obj.LastTorques = Torques;
            obj.LastTorquesx = x;
            if indTorqueMean > 0
                r(indTorqueMean) = Tmean;
                if nargout == 2
                    gradr(indTorqueMean, :) = mean(dTorquesdx, 1);
                end
            end
            if indTorqueStd > 0
                N = numel(obj.OptimizationAngles);
                r(indTorqueStd) = Tstd;
                if nargout == 2
                    gradr(indTorqueStd, :) = 1/Tstd*(1/N*sum(Torques.*dTorquesdx, 1) - Tmean*mean(dTorquesdx, 1));
                end
            end
        end

        function stop = outputFunction(obj, xnorm, optimValues, state)
            stop = false;
            x = xnorm.*(obj.UpperBounds-obj.LowerBounds) + obj.LowerBounds;
            %assert(all(x==obj.LastTorquesx));
            if ~all(x==obj.currentx)
                warning("There might be something wrong here?!")
            end

            obj.Coupling.setRotationAngle(0);
            
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
            obj.OptimizationHistory(pos).Responses(end+1).Name = "Torques";
            obj.OptimizationHistory(pos).Responses(end).Value = obj.LastTorques;

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

            % Display the values of the objectives "Cost" and "TorqueSTD"
            costValue = r(vertcat(obj.OptimizationResponses.Name) == "Cost");
            torqueSTDValue = r(vertcat(obj.OptimizationResponses.Name) == "TorqueSTD");
            torqueMEANValue = r(vertcat(obj.OptimizationResponses.Name) == "TorqueMEAN");
            fprintf('Iteration %d: Cost = %.6f, TorqueSTD = %.6f, TorqueMEAN = %.6f \n', optimValues.iteration, costValue, torqueSTDValue,torqueMEANValue);

            OptimizationHistory = obj.OptimizationHistory;
            save(obj.SaveDirectory+"/OptimizationHistory.mat", 'OptimizationHistory'); 
        end
    
    end
end
