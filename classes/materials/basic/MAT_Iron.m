classdef MAT_Iron < MAT_Magnetic & MAT_Thermal & MAT_Mechanical

    methods (Access = public)
        function obj = MAT_Iron()
            obj.Mur = 500;
            obj.Sigma = 10e6;
            obj.Rho = 7874;
            obj.Cp = 449;
            obj.Kappa = 80;
            obj.PlotColor = [0.7, 0.7, 0.7];
            obj.Nu = 0.3;
            obj.E = 212e9;
        end
    end
end
