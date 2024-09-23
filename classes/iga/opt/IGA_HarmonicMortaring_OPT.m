classdef IGA_HarmonicMortaring_OPT < OPT_Element & IGA_HarmonicMortaring

    properties
        % TBD make list of available responses and check their
        % avalilability when adding responses??!
    end

    methods
        function obj = IGA_HarmonicMortaring_OPT(region1, sides1, region2, sides2)
            obj@IGA_HarmonicMortaring(region1, sides1, region2, sides2)
        end

        function optimizationUpdate(obj, x)
            % do nothing, mortaring is not changed
        end

        function f = f_responses(obj, t, x, y)
            f = zeros(obj.NumberOptimizationResponses, 1);

            a2 = y(obj.Region1.NumberDOF + (1:obj.Region2.NumberDOF), :);
            lambda = y((obj.Region1.NumberDOF+obj.Region2.NumberDOF+1):end, :);

            indTorque = find(vertcat(obj.OptimizationResponses.Name) == "Torque");
            if indTorque > 0
                f(indTorque) = -a2'*obj.CouplingMatrixInit2*obj.RotationMatrixDer*lambda;
            end
        end

        function dfdy = partial_f_partial_y(obj, t, x, y)
            dfdy = sparse(obj.NumberOptimizationResponses, obj.NumberDOF + obj.NumberDOFext);

            a2 = y(obj.Region1.NumberDOF + (1:obj.Region2.NumberDOF), :);
            lambda = y((obj.Region1.NumberDOF+obj.Region2.NumberDOF+1):end, :);

            indTorque = find(vertcat(obj.OptimizationResponses.Name) == "Torque");
            if indTorque > 0
                dfdy(indTorque, :) = -[sparse(obj.Region1.NumberDOF, 1); ...
                    obj.CouplingMatrixInit2*obj.RotationMatrixDer*lambda; ...
                    obj.RotationMatrixDer'*obj.CouplingMatrixInit2'*a2];
            end
        end

        function dfdx = partial_f_partial_x(obj, t, x, y)
            dfdx = sparse(obj.NumberOptimizationResponses, obj.NumberOptimizationParameters);
        end

        function dedx = partial_e_partial_x(obj, t, x, y)
            dedx = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberOptimizationParameters);
        end

        function Lyy = function_Lyy(obj, t, x, y)
            Lyy = sptensor([],[],[obj.NumberDOF, obj.NumberDOF, obj.NumberOptimizationResponses]);

            indTroque = find(vertcat(obj.OptimizationResponses.Name) == "Torque");
            if indTroque > 0
                LyyTorque = zeros(obj.NumberDOF, obj.NumberDOF);
                LyyTorque(obj.a2DOFs, obj.lamdaDOFs) = - obj.CouplingMatrixInit2*obj.RotationMatrixDer;
                LyyTorque(obj.lamdaDOFs, obj.a2DOFs) = - obj.RotationMatrixDer'*obj.CouplingMatrixInit2';
                Lyy(:,:, indTroque) = sptensor(LyyTorque);
            end
        end

        function Lxy = function_Lxy(obj, t, x, y, z)
            Lxy = [];
            % Lxy = sptensor([],[],[obj.NumberOptimizationParameters, obj.NumberDOF, obj.NumberOptimizationResponses]);
        end

        function Lxx = function_Lxx(obj, t, x, y, z)
            Lxx = [];
            % Lxx = sptensor([],[],[obj.NumberOptimizationParameters, obj.NumberOptimizationParameters, obj.NumberOptimizationResponses]);
        end

        function plotFunction(obj, x, optimValues, state)
            % do nothing
        end
    end
end

