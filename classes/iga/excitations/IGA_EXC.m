classdef IGA_EXC < EXC
    properties (Access = public)
        RHS;
    end

    methods
        function obj = IGA_EXC(parentRegion)
            obj@EXC(parentRegion);
        end

        function setValue(obj, func, patches)
            obj.RHS = zeros(obj.ParentRegion.NumberDOF, 1);

            for iPatch = 1:numel(patches)
                patch = patches(iPatch);
                patchDOFs = obj.ParentRegion.Spaces.gnum{patch};
                obj.RHS(patchDOFs) = obj.RHS(patchDOFs) +  op_f_v1(obj.ParentRegion.SpacesEval{patch}, obj.ParentRegion.MeshesEval{patch}, func);
            end
        end

        function f = getEXCvalues(obj, t)
            f = obj.RHS;
        end

    end
end