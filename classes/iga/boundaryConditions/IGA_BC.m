classdef IGA_BC < handle
    properties (Access = public)
        PlotColor;
        PlotAlpha = 1;
        BoundaryNumbers = [];
    end

    methods (Access = public)
        function col = getPlotColor(obj)
            col = obj.PlotColor;
        end

        function alpha = getPlotAlpha(obj)
            alpha = obj.PlotAlpha;
        end
    end
end