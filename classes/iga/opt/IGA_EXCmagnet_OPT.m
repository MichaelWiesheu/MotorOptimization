classdef IGA_EXCmagnet_OPT < IGA_EXCmagnet & OPT_Element
    properties (Access = public)

    end

    methods
        function obj = IGA_EXCmagnet_OPT(parentElement, material, patches)
            obj@IGA_EXCmagnet(parentElement, material, patches);
        end

        function optimizationUpdate(obj, x)
            % if angle of magnet depends on a parameter, update angle
            if isprop(obj.Material, "AngleParameter")
                % Find angle index
                indAngle = vertcat(obj.ParentRegion.OptimizationParameters.Name) == obj.Material.AngleParameter;
                obj.Material.Angle = obj.Material.AngleFunction(x(indAngle));
            end
            obj.updateExcitation();
        end

        function dfdy = partial_f_partial_y(obj, t, x, y)
            dfdy = sparse(obj.ParentRegion.NumberOptimizationResponses, obj.ParentRegion.NumberDOF);
        end

        function dedx = partial_e_partial_x(obj, t, x, y)
            dedx = sparse(obj.ParentRegion.NumberDOF, obj.ParentRegion.NumberOptimizationParameters);
            
            dedx = dedx + op_D_DP_gradv_Br_bot_mp_eval(obj.ParentRegion.Spaces, obj.ParentRegion.SpacesEval, ...
                obj.ParentRegion.SpacesGeo, obj.ParentRegion.SpacesGeoEval, obj.ParentRegion.Meshes, obj.ParentRegion.MeshesEval, ...
                obj.ParentRegion.dCdP, obj.Material, obj.Patches);

            if isprop(obj.Material, "AngleParameter")
                % Find angle index
                indAngle = vertcat(obj.ParentRegion.OptimizationParameters.Name) == obj.Material.AngleParameter;
                angle = obj.Material.AngleFunction(x(indAngle));
                % Numerical evaluation how parameter changes angle
                Step = (obj.ParentRegion.OptimizationParameters(indAngle).Max - obj.ParentRegion.OptimizationParameters(indAngle).Min)*1e-4;
                dAngledP = (obj.Material.AngleFunction(x(indAngle)+Step) - angle)/Step;

                dedx(:, indAngle) = dedx(:, indAngle) + dAngledP.* ...
                    op_d_dalpha_gradv_Br_bot_mp_eval(obj.ParentRegion.Spaces, obj.ParentRegion.SpacesEval, ...
                        obj.ParentRegion.MeshesEval, obj.Material, obj.Patches);
            end            
        end

        function dfdx = partial_f_partial_x(obj, t, x, y)
            % obj.optimizationUpdate(x)
            dfdx = sparse(obj.ParentRegion.NumberOptimizationResponses, obj.ParentRegion.NumberOptimizationParameters);
        end
        
        function f = f_responses(obj, x)
            f = zeros(obj.ParentRegion.NumberOptimizationResponses, 1);
        end

    end
end