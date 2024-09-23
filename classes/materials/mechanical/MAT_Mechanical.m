classdef MAT_Mechanical < MAT
    properties
        IsLinearMECH=true;
        Nu
        E
    end

    methods
        function E = getE(obj)
            E = obj.E;
        end

        function Nu = getNu(obj)
            Nu = obj.Nu;
        end

        function lambda = getLameLambda(obj)
            lambda = obj.getNu()*obj.getE()/((1+obj.getNu())*(1-2*obj.getNu()));
        end

        function mu = getLameMu(obj)
            mu = obj.getE()/(2*(1+obj.getNu()));
        end
    end
end