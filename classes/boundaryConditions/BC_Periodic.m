classdef BC_Periodic < BC
    properties (Access = public)
        LeftDOFs, 
        RightDOFs;
    end

    methods (Access = public)
        function obj = BC_Periodic(parentElement, leftDOFs, rightDOFs)
            leftDOFs = reshape(leftDOFs, 1, []); rightDOFs = reshape(rightDOFs, 1, []);
            assert(isempty(intersect(leftDOFs, rightDOFs)), "Left and right DOFs must not be identical!");
            % Reduce DOFs for possible multiplicity with other Boundary Conditions
            obj@BC(parentElement, union(leftDOFs, rightDOFs));
            obj.LeftDOFs = intersect(obj.DependentDOFs, leftDOFs);
            obj.RightDOFs = intersect(obj.DependentDOFs, rightDOFs);
            % only set LeftDOFs as the dependent DOFs
            obj.DependentDOFs = obj.LeftDOFs;
            % all others are independent
            obj.IndependentDOFs = setdiff((1:obj.ParentElement.NumberDOF), obj.DependentDOFs);
            obj.Name = "Periodic";

            C = sparse(repmat(1:numel(obj.LeftDOFs), 1, 2), [obj.LeftDOFs, obj.RightDOFs], ...
                [ones(1, numel(obj.LeftDOFs)), -ones(1, numel(obj.RightDOFs))], ...
                numel(obj.DependentDOFs), numel(obj.DependentDOFs)+numel(obj.IndependentDOFs));
            obj.CI = C(:, obj.IndependentDOFs);
            obj.CD = C(:, obj.DependentDOFs);
        end

        function [b, db]  = getBCvalues(obj, t)
            b = sparse(numel(obj.DependentDOFs), 1);
            db = sparse(numel(obj.DependentDOFs), 1);
        end
    end
end