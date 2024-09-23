classdef MAT_Aluminum < MAT_Magnetic & MAT_Thermal

    methods (Access = public)
        function obj = MAT_Aluminum()
            obj.Mur = 1;
            obj.Sigma = 35e6;
            obj.Rho = 2700;
            obj.Cp = 896;
            obj.Kappa = 220;
            obj.PlotColor = "#BCC6CC";
            obj.Cost = 2;
        end
    end
end