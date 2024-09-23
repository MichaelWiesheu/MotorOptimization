classdef MAT_Thermal < handle
    properties (Access = protected)
        IsLinearTH=true;
        Cp
        Kappa
    end

    methods (Access = public)
        function Cp = getCp(obj)
            Cp = obj.Cp;
        end
        function cp = setCp(obj, cp)
            obj.Cp = cp;
        end
        
        function kappa = getKappa(obj)
            kappa = obj.Kappa;
        end
        function kappa = setKappa(obj, kappa)
            obj.Kappa = kappa;
        end
    end
end