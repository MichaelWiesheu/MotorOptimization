classdef BC < handle
    properties (Access = public)
        Name;
        ParentElement;
        DependentDOFs;
        IndependentDOFs;
        C
        CI
        CD
    end

    methods (Abstract)
        [b, db]  = getBCvalues(obj, t)
    end

    methods (Access = public)
        function obj = BC(parentElement, dependentDOFs)
            obj.ParentElement = parentElement;
            obj.DependentDOFs = reshape(dependentDOFs, 1, []);
            % Check all other boundary conditions for consistency
            for iBC = 1:numel(obj.ParentElement.BoundaryConditions)
                checkBC = obj.ParentElement.BoundaryConditions(iBC).BoundaryCondition;
                occupiedDOFs = intersect(obj.DependentDOFs, checkBC.DependentDOFs, 'stable');
                if ~isempty(occupiedDOFs)
                    obj.DependentDOFs = setdiff(obj.DependentDOFs, occupiedDOFs, 'stable');
                    disp(['Note: DOF(s) ' num2str(reshape(occupiedDOFs, 1, [])) ' were ignored in setting the boundary condition since they are alreay in use.'])
                end
            end
            obj.IndependentDOFs = setdiff((1:obj.ParentElement.NumberDOF), obj.DependentDOFs);
        end

        function [dofsD, dofsI]  = getBCindices(obj)
            dofsD = obj.DependentDOFs;
            dofsI = obj.IndependentDOFs;
        end

        function [G, H] = getBCmatrices(obj)
            G = -inv(obj.CD)*obj.CI;
            H = inv(obj.CD);
        end
    end
end