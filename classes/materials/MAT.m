classdef MAT < handle

    properties
        Rho=0;
        Cost=0;
        PlotColor= [1, 1, 1];
        PlotAlpha = 1;
    end

    methods (Access = public)
        
        function rho = getRho(obj)
            rho = obj.Rho;
        end
        function rho = setRho(obj, rho)
            obj.Rho = rho;
        end

        function cost = getCost(obj)
            cost = obj.Cost;
        end
        function setCost(obj, cost)
            obj.Cost = cost;
        end

        function col = getPlotColor(obj)
            col = obj.PlotColor;
        end

        function alpha = getPlotAlpha(obj)
            alpha = obj.PlotAlpha;
        end
    end
end