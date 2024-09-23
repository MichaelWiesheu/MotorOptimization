classdef IGA_BCantiPeriodic < IGA_BC & BC_AntiPeriodic

    methods (Access = public)
        function obj = IGA_BCantiPeriodic(parentElement, leftSides, rightSides)
            
            leftDOFs = parentElement.getBoundaryDOFs(leftSides);
            rightDOFs = parentElement.getBoundaryDOFs(rightSides);

            obj@BC_AntiPeriodic(parentElement, leftDOFs, rightDOFs);

            obj.BoundaryNumbers = union(leftSides, rightSides);
            obj.PlotColor = TUDa_getColor("5a");
        end
    end
end