classdef MAT_Copper < MAT_Magnetic & MAT_Thermal & MAT_Mechanical

    methods (Access = public)
        function obj = MAT_Copper()
            obj.Mur = 1;
            obj.Sigma = 5.96e7;
            obj.Rho = 8960;
            obj.Cp = 382;
            obj.Kappa = 401;
            obj.PlotColor = "#E9503E";
            obj.Cost = 10;
            % https://www.mit.edu/~6.777/matprops/copper.htm
            obj.Nu = 0.34; 
            obj.E = 13e7;
        end
    end
end