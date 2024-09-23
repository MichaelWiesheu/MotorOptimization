classdef IGA_BCdirichlet < IGA_BC & BC_Dirichlet
    properties (Access = public)
        Spaces, SpacesEval;
        Meshes, MeshesEval;
        const_f;
    end

    methods (Access = public)
        function obj = IGA_BCdirichlet(parentElement, boundaryNumbers, value, dvalue)
            % for transient problems, derivative needs to be given: 
            % value = @(t) @(x,y) function of size x, depending of t, x, y
            % dvalue = @(t) @(x,y) function of size x, depending of t, x, y
            % for static problems:
            % value = @(x,y) function of size x
            % value = const
            if ~exist("value", "var")
                value = 0;
            end
            if ~exist("dvalue", "var")
                dvalue = 0;
            end
            dDOFs = parentElement.getBoundaryDOFs(boundaryNumbers);
            
            obj@BC_Dirichlet(parentElement, dDOFs, value, dvalue);
            obj.BoundaryNumbers = boundaryNumbers;
            obj.PlotColor = TUDa_getColor("2b");
            
            % Generate Spaces and Meshes
            for iBnd = 1:numel(boundaryNumbers)
                patch = parentElement.Boundaries(boundaryNumbers(iBnd)).patches;
                side = parentElement.Boundaries(boundaryNumbers(iBnd)).faces;
                obj.MeshesEval{iBnd} = msh_eval_boundary_side(parentElement.Meshes.msh_patch{patch}, side);
                msh_side_int = msh_boundary_side_from_interior(parentElement.Meshes.msh_patch{patch}, side);
                obj.Spaces{iBnd} = parentElement.Spaces.sp_patch{patch}.constructor(msh_side_int);
                obj.SpacesEval{iBnd} = sp_precompute(obj.Spaces{iBnd}, msh_side_int, 'value', true, 'gradient', true);
            end

            % Generate IGA projection matrix:
            M = sparse(parentElement.NumberDOF, parentElement.NumberDOF);
            obj.const_f = sparse(parentElement.NumberDOF, 1);       % precalculate intetral for spatially constant values
            for iBnd = 1:numel(boundaryNumbers)
                patch = parentElement.Boundaries(boundaryNumbers(iBnd)).patches;
                patchDOFs = parentElement.Spaces.gnum{patch};
                M(patchDOFs, patchDOFs) = M(patchDOFs, patchDOFs) + ...
                    + 1/parentElement.Length*op_u_v (obj.SpacesEval{iBnd}, obj.SpacesEval{iBnd}, obj.MeshesEval{iBnd}, 1);
                obj.const_f(patchDOFs) = obj.const_f(patchDOFs) + op_f_v1(obj.SpacesEval{iBnd}, obj.MeshesEval{iBnd});
            end
            
            obj.CD = M(obj.DependentDOFs, obj.DependentDOFs);
            obj.CI = sparse(numel(obj.DependentDOFs), numel(obj.IndependentDOFs));
        end

        function [b, db] = getBCvalues(obj, t)
            if isa(obj.ValueDerivative, "function_handle")
               % Transient problem
               func = obj.Value(t);
               funcDer = obj.ValueDerivative(t);
               if isa(func, "function_handle")
                    % Spatially changing value:
                    b = sparse(obj.ParentElement.NumberDOF, 1);
                    for iBnd = 1:numel(obj.BoundaryNumbers)
                        patch = obj.ParentElement.Boundaries(obj.BoundaryNumbers(iBnd)).patches;
                        patchDOFs = obj.ParentElement.Spaces.gnum{patch};
                        b(patchDOFs) = b(patchDOFs) + op_f_v1(obj.SpacesEval{iBnd}, obj.MeshesEval{iBnd}, func);
                    end
                    b = b(obj.DependentDOFs);
                    if nargout == 2
                        db = sparse(obj.ParentElement.NumberDOF, 1);
                        for iBnd = 1:numel(obj.BoundaryNumbers)
                            patch = obj.ParentElement.Boundaries(obj.BoundaryNumbers(iBnd)).patches;
                            patchDOFs = obj.ParentElement.Spaces.gnum{patch};
                            db(patchDOFs) = db(patchDOFs) + op_f_v1(obj.SpacesEval{iBnd}, obj.MeshesEval{iBnd}, funcDer);
                        end
                        db = db(obj.DependentDOFs);
                    end
                else
                    % Spatially constant value:
                    b = func*obj.const_f(obj.DependentDOFs);
                    if nargout == 2
                        db = funcDer*obj.const_f(obj.DependentDOFs);
                    end
               end
            else
                % Static problem                
                if isa(obj.Value, "function_handle")
                    % Spatially changing value:
                    b = sparse(obj.ParentElement.NumberDOF, 1);
                    for iBnd = 1:numel(obj.BoundaryNumbers)
                        patch = obj.ParentElement.Boundaries(obj.BoundaryNumbers(iBnd)).patches;
                        patchDOFs = obj.ParentElement.Spaces.gnum{patch};
                        b(patchDOFs) = b(patchDOFs) + op_f_v1(obj.SpacesEval{iBnd}, obj.MeshesEval{iBnd}, obj.Value);
                    end
                    b = b(obj.DependentDOFs);
                    if nargout == 2
                        db = obj.ValueDerivative*obj.const_f(obj.DependentDOFs);
                    end
                else
                    % Spatially constant value:
                    b = obj.Value*obj.const_f(obj.DependentDOFs);
                    if nargout == 2
                        db = obj.ValueDerivative*obj.const_f(obj.DependentDOFs);
                    end
                end
            end

            % if isa(obj.Value, "function_handle")
            %     b = ones(numel(obj.DependentDOFs),1)*obj.Value(t);
            % else
            %     b = ones(numel(obj.DependentDOFs),1)*obj.Value;
            % end
            % 
            % if isa(obj.ValueDerivative, "function_handle")
            %     db = ones(numel(obj.DependentDOFs),1)*obj.ValueDerivative(t);
            % else
            %     db = ones(numel(obj.DependentDOFs),1)*obj.ValueDerivative;
            % end
        end
    end
end