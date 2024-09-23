classdef EXC < handle
    properties (Access = public)
        Name;
        ParentRegion;
        IsConstant=false;
    end

    methods (Abstract)
        [f, df]  = getEXCvalues(obj, t)
    end

    methods (Access = public)
        function obj = EXC(parentRegion)
            obj.ParentRegion = parentRegion;
        end
    end
end