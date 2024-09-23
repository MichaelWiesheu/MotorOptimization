classdef IGA_RegionMagnetic_OPT <  IGA_RegionMagnetic & IGA_Region_OPT

    properties

    end

    methods
        function obj = IGA_RegionMagnetic_OPT()

        end

        function setMaterial(obj, material, patches)
            setMaterial@IGA_Region(obj, material, patches);
            if isprop(material, 'Br') && isprop(material, 'Angle')
                MagnetExcitation = IGA_EXCmagnet_OPT(obj, material, patches);
                obj.addExcitation(MagnetExcitation);
            end
        end

        function dedx = partial_e_partial_x(obj, t, x, y)
            dedx = sparse(obj.NumberDOF, obj.NumberOptimizationParameters);
            % calculate only, if geometry changes
            if ~isempty(obj.GeometryFunction)
                for iMat = 1:numel(obj.Materials)
                    patches = obj.Materials(iMat).Patches;
                    mat = obj.Materials(iMat).Material;
                    if obj.Materials(iMat).Material.IsLinearMAG == true
                        nu = obj.Materials(iMat).Material.getNuLinear();
                        % Length scaling in multiplication
                        dedx = dedx - nu*op_D_DP_gradu_gradv_times_u_mp_eval(obj.Spaces, obj.SpacesEval, obj.Spaces, obj.SpacesEval, ...
                            obj.SpacesGeo, obj.SpacesGeoEval, obj.Meshes, obj.MeshesEval, obj.dCdP, y/obj.Length, patches);
                    elseif obj.Materials(iMat).Material.IsLinearMAG == false
                        % Length scaling in multiplication
                        dedx = dedx - op_D_DP_gradu_nu_gradv_times_u_mp_eval(obj.Spaces, obj.SpacesEval, obj.Spaces, obj.SpacesEval, ...
                            obj.SpacesGeo, obj.SpacesGeoEval, obj.Meshes, obj.MeshesEval, obj.dCdP, y/obj.Length, mat, patches);
                    end
                end
            end

            for iEXC = 1:numel(obj.Excitations)
                if isa(obj.Excitations(iEXC).Excitation, "OPT_Element")
                    dedx = dedx + obj.Excitations(iEXC).Excitation.partial_e_partial_x(t, x, y);
                end
            end
        end

    end
end

