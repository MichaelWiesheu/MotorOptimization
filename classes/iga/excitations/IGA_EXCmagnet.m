classdef IGA_EXCmagnet < EXC
    properties (Access = public)
        RHS;
        Material;
        Patches;
    end

    methods
        function obj = IGA_EXCmagnet(parentRegion, material, patches)
            obj@EXC(parentRegion);
            obj.IsConstant = true;
            obj.Material = material;
            obj.Patches = patches;
            obj.updateExcitation();
        end

        function f = getEXCvalues(obj, t)
            f = obj.RHS;
        end

        function updateExcitation(obj)
            obj.RHS = op_gradv_Br_bot_mp(obj.ParentRegion.Spaces, obj.ParentRegion.Meshes, obj.Material, obj.Patches);
        end

        function plotExcitation(obj)
            % Plotting functions here for arrows and so!
        end
    end
end