classdef IGA_BCperiodicNormal < IGA_BC & BC_PeriodicNormal

    methods (Access = public)
        function obj = IGA_BCperiodicNormal(parentElement, leftSides, rightSides, nLeft, nRight)
            
            [leftDOFsX, leftDOFsY] = parentElement.getBoundaryDOFs(leftSides);
            [rightDOFsX, rightDOFsY]  = parentElement.getBoundaryDOFs(rightSides);

            obj@BC_PeriodicNormal(parentElement, leftDOFsX, leftDOFsY, rightDOFsX, rightDOFsY, nLeft, nRight);

            obj.BoundaryNumbers = union(leftSides, rightSides);
            obj.PlotColor = TUDa_getColor("11c");
        end
    end
end

