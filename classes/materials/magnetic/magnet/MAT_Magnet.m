classdef MAT_Magnet < MAT_Magnetic & MAT_Thermal & MAT_Mechanical
    properties (Access = public)
        Br;
        Angle;
    end
    methods (Access = public)
        function obj = MAT_Magnet(angle)
            obj.Angle = angle;
            % Standard values
            obj.Cost = 50;
            obj.Mur = 1.05;
            obj.Sigma = 1/(1.4e-6);
            obj.Rho = 7500;
            obj.PlotColor = "#00715E";
            obj.Nu = 0.3;
            obj.E = 180e9;
         end
    end
end