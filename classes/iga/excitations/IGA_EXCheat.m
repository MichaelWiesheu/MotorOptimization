classdef IGA_EXCheat < EXC
    properties (Access = public)
        RHS;
        Flux;
        Patches;
    end

    methods
        function obj = IGA_EXCheat(parentRegion, patches, fluxFunction)
            obj@EXC(parentRegion);
            obj.Flux = fluxFunction;
            if isa(fluxFunction, "function_handle")
                obj.IsConstant = false;
            else
                obj.IsConstant = true;
            end

            obj.Patches = patches;
            obj.updateExcitation()
        end

        function updateExcitation(obj)
            obj.RHS = op_f_v_mp_eval (obj.ParentRegion.Spaces, obj.ParentRegion.SpacesEval, obj.ParentRegion.MeshesEval, obj.Patches);
        end

        function f = getEXCvalues(obj, t)
            if obj.IsConstant
                f = obj.Flux*obj.RHS;
            else
                f = obj.Flux(t)*obj.RHS;
            end
        end

    end
end