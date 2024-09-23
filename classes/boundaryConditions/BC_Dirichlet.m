classdef BC_Dirichlet < BC
    properties (Access = public)
        Value
        ValueDerivative
    end

    methods (Access = public)
        function obj = BC_Dirichlet(parentElement, dDOFs, value, dvalue)
            obj@BC(parentElement, dDOFs);
            obj.Name = "Dirichlet";
            if ~exist("value", "var")
                value = 0;
            end
            if ~exist("dvalue", "var")
                dvalue = 0;
            end
            obj.ValueDerivative = dvalue;
            obj.Value = value;

            obj.CI = sparse(numel(obj.DependentDOFs), numel(obj.IndependentDOFs));
            obj.CD = sparse(eye(numel(obj.DependentDOFs)));
        end

        function [b, db] = getBCvalues(obj, t)
            if isa(obj.Value, "function_handle")
                b = ones(numel(obj.DependentDOFs),1)*obj.Value(t);
            else
                b = ones(numel(obj.DependentDOFs),1)*obj.Value;
            end

            if isa(obj.ValueDerivative, "function_handle")
                db = ones(numel(obj.DependentDOFs),1)*obj.ValueDerivative(t);
            else
                db = ones(numel(obj.DependentDOFs),1)*obj.ValueDerivative;
            end
        end
    end
end