classdef MAT_NdFeB_OPT < MAT_Magnet
    properties (Access = public)
        AngleFunction;
        AngleParameter;
    end

    methods (Access = public)
        function obj = MAT_NdFeB_OPT(Br, angle, angleFunction, angleParameter)
            obj@MAT_Magnet(angle);
            obj.Mur = 1.05;
            obj.Br = Br;
            obj.AngleFunction = angleFunction;  % Function for the magnet orientation depending on parameter
            obj.AngleParameter = angleParameter;  % Parameter name for the orientation parameter
         end
    end
end