classdef MAT_Air < MAT_Magnetic & MAT_Thermal & MAT_Mechanical

    methods (Access = public)
        function obj = MAT_Air()
            %https://theengineeringmindset.com/properties-of-air-at-atmospheric-pressure/
            obj.Mur = 1;
            obj.Sigma = 0;
            obj.Rho = 1.2047;
            obj.Cp = 1006.1;
            obj.Kappa = 0.025596;
            obj.PlotColor = "#009CDA";

            obj.Nu = 0.3;
            obj.E = 1e5;
            obj.Cost = 0;
        end
    end
end