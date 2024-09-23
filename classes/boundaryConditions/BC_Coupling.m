classdef BoundaryCoupling < BoundaryCondition
    properties (Access = public)

    end

    methods (Access = public)
        function obj = BoundaryCoupling(parentElement, dDOFs)
            obj@BoundaryCondition(parentElement, dDOFs);
            obj.Name = "Coupling";
            obj.PlotColor= TUDa_getColor("9b");
        end

        function [dofsD, dofsI, G, H, b]  = getBoundaryMatrices(obj, t)
            dofsD = obj.DependentDOFs;
            dofsI = obj.IndependentDOFs;
            CI = sparse(0, numel(obj.IndependentDOFs));
            CD = sparse(0, 0);
            G = -CD*CI; %inv(CD)
            H = CD; % inv(CD)
            b = sparse(0, 1);
        end
    end
end