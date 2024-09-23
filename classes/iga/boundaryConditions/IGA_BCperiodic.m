classdef IGA_BCperiodic < IGA_BC & BC_Periodic

    methods (Access = public)
        function obj = IGA_BCperiodic(parentElement, leftSides, rightSides)
            
            leftDOFs = parentElement.getBoundaryDOFs(leftSides);
            rightDOFs = parentElement.getBoundaryDOFs(rightSides);

            obj@BC_Periodic(parentElement, leftDOFs, rightDOFs);

            obj.BoundaryNumbers = union(leftSides, rightSides);
            obj.PlotColor = TUDa_getColor("11a");
        end
    end
end