classdef MAT_NdFeB_Br14 < MAT_Magnet
    
    methods (Access = public)
        function obj = MAT_NdFeB_Br14(angle)
            obj@MAT_Magnet(angle);
            obj.Mur = 1.05;
            obj.Br = 1.4;
         end
    end
end