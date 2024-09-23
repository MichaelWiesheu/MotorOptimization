classdef IGA_HarmonicMortaring < DAE_Element
    properties (Access = public)
        Boundaries1, Boundaries2;
        Region1, Region2;
        HarmonicsSin, HarmonicsCos, HarmonicsAll;
        %         NumberDOF;
        CouplingIndices;
        RotationAngle = 0;
        Meshes1, Spaces1, MeshesEval1, SpacesEval1;
        Meshes2, Spaces2, MeshesEval2, SpacesEval2;
        CouplingMatrixFull1, CouplingMatrixInit1, CouplingMatrix1;
        CouplingMatrixFull2, CouplingMatrixInit2, CouplingMatrix2;
        %         RotationMatrix1, RotationMatrix2;
        %         RotationAngle1, RotationAngle2;
        RotationMatrix, RotationMatrixDer, RotationMatrixDer2;
        MaterialParameter = 1/(4*pi*1e-7);
        a1DOFs;
        a2DOFs;
        lamdaDOFs;
    end

    methods (Access = public)
        function obj = IGA_HarmonicMortaring(region1, sides1, region2, sides2, materialParameter)
            if exist("materialParameter", "var")
                obj.MaterialParameter = materialParameter;
            end
            obj.addExternalElement(region1)
            obj.addExternalElement(region2)

            obj.Region1 = region1;
            obj.Region2 = region2;
            obj.Boundaries1 = sides1;
            obj.Boundaries2 = sides2;
        end

        function setHarmonics(obj, sinValues, cosValues)
            obj.HarmonicsSin = reshape(unique(sinValues), 1, []);
            obj.HarmonicsCos = reshape(unique(cosValues), 1, []);
            obj.HarmonicsAll = 0:max([obj.HarmonicsSin, obj.HarmonicsCos]);
            obj.NumberDOF = numel(obj.HarmonicsSin) + numel(obj.HarmonicsCos);
            obj.a1DOFs = 1:obj.Region1.NumberDOF;
            obj.a2DOFs = obj.Region1.NumberDOF+(1:obj.Region2.NumberDOF);
            obj.lamdaDOFs = obj.NumberDOFext+(1:obj.NumberDOF);
            obj.CouplingIndices = union(obj.HarmonicsSin*2 + 1, obj.HarmonicsCos*2 + 2);
            obj.generateSpacesAndMeshes();
            obj.generateCouplingMatrices();
        end

        function generateSpacesAndMeshes(obj)
            for iBnd1 = 1:numel(obj.Boundaries1)
                patch = obj.Region1.Boundaries(obj.Boundaries1(iBnd1)).patches;
                side = obj.Region1.Boundaries(obj.Boundaries1(iBnd1)).faces;
                obj.MeshesEval1{iBnd1} = msh_eval_boundary_side(obj.Region1.Meshes.msh_patch{patch}, side);
                msh_side_int = msh_boundary_side_from_interior(obj.Region1.Meshes.msh_patch{patch}, side);
                obj.Spaces1{iBnd1} = obj.Region1.Spaces.sp_patch{patch}.constructor(msh_side_int);
                obj.SpacesEval1{iBnd1} = sp_precompute(obj.Spaces1{iBnd1}, msh_side_int, 'value', true, 'gradient', true);
            end
            % region2
            for iBnd2 = 1:numel(obj.Boundaries2)
                patch = obj.Region2.Boundaries(obj.Boundaries2(iBnd2)).patches;
                side = obj.Region2.Boundaries(obj.Boundaries2(iBnd2)).faces;
                obj.MeshesEval2{iBnd2} = msh_eval_boundary_side(obj.Region2.Meshes.msh_patch{patch}, side);
                msh_side_int = msh_boundary_side_from_interior(obj.Region2.Meshes.msh_patch{patch}, side);
                obj.Spaces2{iBnd2} = obj.Region2.Spaces.sp_patch{patch}.constructor(msh_side_int);
                obj.SpacesEval2{iBnd2} = sp_precompute(obj.Spaces2{iBnd2}, msh_side_int, 'value', true, 'gradient', true);
            end
        end

        function generateCouplingMatrices(obj)
            progressbar('Generating coupling matrices: ')
            % define array functions for all harmonics
            SinPhi =  @(x, y) sin(reshape(obj.HarmonicsAll, [], 1).* atan2 (reshape(y, [1, size(y)]), reshape(x, [1, size(x)])));
            CosPhi =  @(x, y) cos(reshape(obj.HarmonicsAll, [], 1).* atan2 (reshape(y, [1, size(y)]), reshape(x, [1, size(x)])));

            % region1
            B1 = sparse (obj.Region1.NumberDOF, 2*numel(obj.HarmonicsAll));
            for iBnd = 1:numel(obj.Boundaries1)
                patchDOFs = obj.Region1.Spaces.gnum{obj.Region1.Boundaries(obj.Boundaries1(iBnd)).patches};
                B1(patchDOFs, 1:2:end) = B1(patchDOFs, 1:2:end) + op_fs_v (obj.SpacesEval1{iBnd}, obj.MeshesEval1{iBnd}, SinPhi);
                B1(patchDOFs, 2:2:end) = B1(patchDOFs, 2:2:end) + op_fs_v (obj.SpacesEval1{iBnd}, obj.MeshesEval1{iBnd}, CosPhi);
                progressbar(iBnd/(numel(obj.Boundaries1)+numel(obj.Boundaries2)));
            end

            % region2
            B2 = sparse (obj.Region2.NumberDOF, 2*numel(obj.HarmonicsAll));
            for iBnd = 1:numel(obj.Boundaries2)
                patchDOFs = obj.Region2.Spaces.gnum{obj.Region2.Boundaries(obj.Boundaries2(iBnd)).patches};
                B2(patchDOFs, 1:2:end) = B2(patchDOFs, 1:2:end) + op_fs_v (obj.SpacesEval2{iBnd}, obj.MeshesEval2{iBnd}, SinPhi);
                B2(patchDOFs, 2:2:end) = B2(patchDOFs, 2:2:end) + op_fs_v (obj.SpacesEval2{iBnd}, obj.MeshesEval2{iBnd}, CosPhi);
                progressbar((numel(obj.Boundaries1)+iBnd)/(numel(obj.Boundaries1)+numel(obj.Boundaries2)));
            end

            B1 = B1 * obj.MaterialParameter;
            B2 = B2 * obj.MaterialParameter;

            obj.CouplingMatrixFull1 = B1;
            obj.CouplingMatrixFull2 = B2;

            obj.CouplingMatrixInit1 = obj.CouplingMatrixFull1(:, obj.CouplingIndices);
            obj.CouplingMatrixInit2 = obj.CouplingMatrixFull2(:, obj.CouplingIndices);

            obj.CouplingMatrix1 = obj.CouplingMatrixInit1;
            obj.CouplingMatrix2 = obj.CouplingMatrixInit2;

            obj.setRotationAngle(obj.RotationAngle);
            progressbar('finished');
        end

        function setRotationAngle(obj, angleRad)
            obj.RotationAngle = angleRad;

            nHarm = numel(obj.HarmonicsAll);
            rows = [1:2:nHarm*2, 1:2:nHarm*2, 2:2:nHarm*2, 2:2:nHarm*2];
            cols = [1:2:nHarm*2, 2:2:nHarm*2, 1:2:nHarm*2, 2:2:nHarm*2];
            vals = [cos(obj.HarmonicsAll.*obj.RotationAngle), sin(obj.HarmonicsAll.*obj.RotationAngle), ...
                   -sin(obj.HarmonicsAll.*obj.RotationAngle), cos(obj.HarmonicsAll.*obj.RotationAngle)];
            R = sparse(rows, cols, vals, 2*nHarm, 2*nHarm);

            valsDer = [-obj.HarmonicsAll.*sin(obj.HarmonicsAll.*obj.RotationAngle), obj.HarmonicsAll.*cos(obj.HarmonicsAll.*obj.RotationAngle), ...
                       -obj.HarmonicsAll.*cos(obj.HarmonicsAll.*obj.RotationAngle), -obj.HarmonicsAll.*sin(obj.HarmonicsAll.*obj.RotationAngle)];
            Rder = sparse(rows, cols, valsDer, 2*nHarm, 2*nHarm);

            valsDer2 = [-obj.HarmonicsAll.^2.*cos(obj.HarmonicsAll.*obj.RotationAngle), -obj.HarmonicsAll.^2.*sin(obj.HarmonicsAll.*obj.RotationAngle), ...
                         obj.HarmonicsAll.^2.*sin(obj.HarmonicsAll.*obj.RotationAngle), -obj.HarmonicsAll.^2.*cos(obj.HarmonicsAll.*obj.RotationAngle)];
            Rder2 = sparse(rows, cols, valsDer2, 2*nHarm, 2*nHarm);

            obj.RotationMatrix = R(obj.CouplingIndices, obj.CouplingIndices);
            obj.RotationMatrixDer = Rder(obj.CouplingIndices, obj.CouplingIndices);
            obj.RotationMatrixDer2 = Rder2(obj.CouplingIndices, obj.CouplingIndices);

            Bst = obj.CouplingMatrixFull2*R;
            obj.CouplingMatrix2 = Bst(:, obj.CouplingIndices);

            obj.Region1.updatePlotPoints(angleRad);
        end

        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            if nargout == 2
                Mlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            % Nonlinear because matrix can change with rotation angle
            Knonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            Knonlin(obj.a2DOFs, obj.lamdaDOFs) = -obj.CouplingMatrix2;
            Knonlin(obj.lamdaDOFs, obj.a2DOFs) = -obj.CouplingMatrix2';
            if nargout == 2
                Klin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
                Klin(obj.a1DOFs, obj.lamdaDOFs) = obj.CouplingMatrix1;
                Klin(obj.lamdaDOFs, obj.a1DOFs) = obj.CouplingMatrix1';
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(obj.NumberDOFext+obj.NumberDOF, 1);
            if nargout == 2
                Flin = sparse(obj.NumberDOFext+obj.NumberDOF, 1);
            end
        end

        function J = JacobianMatrix(obj, t, y)
            J = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
        end

        function postprocess(obj, ts, ys)
            a1 = ys(1:obj.Region1.NumberDOF, :);
            a2 = ys(obj.Region1.NumberDOF + (1:obj.Region2.NumberDOF), :);
            obj.Solution.Values = ys((obj.Region1.NumberDOF+obj.Region2.NumberDOF+1):end, :);
            % NumberDOFext
            if isscalar(ts)
                obj.Solution.Torque = -a2'*obj.CouplingMatrixInit2*obj.RotationMatrixDer*obj.Solution.Values;
            end 
        end

        function T = calcTorqueBrBtRt(obj)
            T = 0;
            a = obj.Region1.MagneticPotential;
            for iBnd = 1:numel(obj.Boundaries1)
                patchDOFs = obj.Region1.Spaces.gnum{obj.Region1.Boundaries(obj.Boundaries1(iBnd)).patches};
                a_loc = a(patchDOFs);
                T = T + obj.MaterialParameter * obj.Region1.Length * ...
                op_Br_Bt_dGamma(obj.SpacesEval1{iBnd}, obj.MeshesEval1{iBnd}, a_loc);
            end
        end

        function T = calcTorqueBrBtSt(obj)
            T = 0;
            a = obj.Region2.MagneticPotential;
            for iBnd = 1:numel(obj.Boundaries2)
                patchDOFs = obj.Region2.Spaces.gnum{obj.Region2.Boundaries(obj.Boundaries2(iBnd)).patches};
                a_loc = a(patchDOFs);
                T = T + obj.MaterialParameter * obj.Region2.Length * ...
                op_Br_Bt_dGamma(obj.SpacesEval2{iBnd}, obj.MeshesEval2{iBnd}, a_loc);
            end
        end

        function [q1, q2] = calcHeatFlux(obj)
            q1 = 0;
            T = obj.Region1.Temperature;
            for iBnd = 1:numel(obj.Boundaries1)
                patchDOFs = obj.Region1.Spaces.gnum{obj.Region1.Boundaries(obj.Boundaries1(iBnd)).patches};
                t_loc = T(patchDOFs);
                q1 = q1 + obj.MaterialParameter * ...
                op_gradTn_dGamma(obj.SpacesEval1{iBnd}, obj.MeshesEval1{iBnd}, t_loc);
            end
            q2 = 0;
            T = obj.Region2.Temperature;
            for iBnd = 1:numel(obj.Boundaries2)
                patchDOFs = obj.Region2.Spaces.gnum{obj.Region2.Boundaries(obj.Boundaries2(iBnd)).patches};
                t_loc = T(patchDOFs);
                q2 = q2 + obj.MaterialParameter * ...
                op_gradTn_dGamma(obj.SpacesEval2{iBnd}, obj.MeshesEval2{iBnd}, t_loc);
            end
        end

        function q = calcHeatFluxMortar(obj)
            q = sum(obj.CouplingMatrix1, 1)*obj.Solution.Values;
            q = ones(1, size(obj.CouplingMatrix1, 1))*obj.CouplingMatrix1*obj.Solution.Values
        end

        function plotCouplingValues(obj)
            figure;
            sinValues = obj.Solution.Values((mod(obj.CouplingIndices, 2) == 1));
            cosValues = obj.Solution.Values((mod(obj.CouplingIndices, 2) == 0));
            xlims = [min(obj.HarmonicsAll)-1, max(obj.HarmonicsAll)+1];
            subplot(2,1,1)
            bar(obj.HarmonicsSin, sinValues);
            xticks(obj.HarmonicsSin);
            xlim(xlims);
            title("Sin Values");
            subplot(2,1,2)
            bar(obj.HarmonicsCos, cosValues);
            xticks(obj.HarmonicsCos);
            xlim(xlims);
            title("Cos Values");
        end
    end
end