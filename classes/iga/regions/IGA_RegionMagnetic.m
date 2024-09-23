classdef IGA_RegionMagnetic < IGA_Region

    properties
        NewtonResidual = 1e-6;
        MagneticPotential;
        FigureRemanence;
        FigureBRes, FigureBResBMag, FigureBResBoundaries, FigureBResInterfaces, FigureBResQuiver
        FigurePotLines;
        FigurePotLinesGeometry
        FigurePotLinesBoundaries
        FigurePotLinesInterfaces
        FigurePotLinesPotentialLines

        PlotValuesNeedUpdate = true;
    end

    methods
        function obj = IGA_RegionMagnetic()
            obj@IGA_Region();
        end

        function generateSpaces(obj)
            progressbar([obj.Name ': Generating function spaces: '])

            % pre-evaluated spaces
            obj.SpacesEval  = cell (1, obj.NumberPatches);

            % temporary cell arrays
            spaces  = cell (1, obj.NumberPatches);
            for iPatch = 1:obj.NumberPatches
                [knots, ~] = kntrefine (obj.Geometry(iPatch).nurbs.knots, obj.SubdivisionsPatches-1, obj.FormFunctionDegree, obj.FormFunctionDegree - 1);
                spaces{iPatch}  = sp_bspline (knots, obj.FormFunctionDegree, obj.Meshes.msh_patch{iPatch});
                obj.SpacesEval{iPatch}  = sp_precompute(spaces{iPatch}, obj.MeshesEval{iPatch}, 'value', true, 'gradient', true);
                progressbar(iPatch/obj.NumberPatches);
            end
            obj.Spaces  = sp_multipatch(spaces, obj.Meshes, obj.Interfaces, obj.BoundaryInterfaces);

            % init values
            obj.NumberDOF = obj.Spaces.ndof;
            obj.MagneticPotential = zeros(obj.NumberDOF, 1);
            obj.PlotValuesNeedUpdate = true;

            progressbar('finished');
        end

        function setMaterial(obj, material, patches)
            setMaterial@IGA_Region(obj, material, patches);
            if isprop(material, 'Br') && isprop(material, 'Angle')
                MagnetExcitation = IGA_EXCmagnet(obj, material, patches);
                obj.addExcitation(MagnetExcitation);
            end
        end

        function plotGeometry(obj)
            plotGeometry@IGA_Region(obj);
            obj.plotRemanence();
        end

        function plotRemanence(obj)
            xvals = cell(obj.NumberPatches, 1);
            yvals = cell(obj.NumberPatches, 1);
            Bx = cell(obj.NumberPatches, 1);
            By = cell(obj.NumberPatches, 1);
            for iMat = 1:numel(obj.Materials)
                if isprop(obj.Materials(iMat).Material, 'Br')
                    for patch = obj.Materials(iMat).Patches
                        xvals{patch} = reshape(obj.Geometry(patch).PlotPointsQuiver{1}, [], 1);
                        yvals{patch} = reshape(obj.Geometry(patch).PlotPointsQuiver{2}, [], 1);
                        Bx{patch} = reshape(ones(size(obj.Geometry(patch).PlotPointsQuiver{1})), [], 1)*cos(obj.Materials(iMat).Material.Angle);
                        By{patch} = reshape(ones(size(obj.Geometry(patch).PlotPointsQuiver{2})), [], 1)*sin(obj.Materials(iMat).Material.Angle);
                    end
                end
            end
            xvals = cell2mat(xvals);
            yvals = cell2mat(yvals);
            Bx = cell2mat(Bx);
            By = cell2mat(By);
            % Rotation
            Bxtemp = Bx*cos(obj.RotationAngle) - By*sin(obj.RotationAngle);
            Bytemp = By*cos(obj.RotationAngle) + Bx*sin(obj.RotationAngle);
            Bx = Bxtemp;
            By = Bytemp;
            if ishandle(obj.FigureRemanence)
                set(obj.FigureRemanence, 'XData', xvals, 'YData', yvals, 'UData', Bx, 'VData', By)
            else
                obj.FigureRemanence = quiver(xvals, yvals, Bx, By, 2e-1, 'color', 'k', 'LineWidth', 1.5);
            end
        end

        function calculatePlotPoints(obj, nsubs, patches)
            if nargin <= 1
                nsubs = [1, 1];
            end
            if nargin <= 2
                patches = 1:obj.NumberPatches;
            end
            calculatePlotPoints@IGA_Region(obj, nsubs, patches);
            obj.calculatePlotPointsQuiver();
        end

        function calculatePlotPointsQuiver(obj, nQuivers)
            if ~exist("nQuivers", "var")
                nQuivers = 300;
            end
            % lengths of patch sides
            for iPatch = 1:obj.NumberPatches
                l1 = sum(msh_eval_boundary_side(obj.Meshes.msh_patch{iPatch}, 1).element_size);
                l2 = sum(msh_eval_boundary_side(obj.Meshes.msh_patch{iPatch}, 2).element_size);
                l3 = sum(msh_eval_boundary_side(obj.Meshes.msh_patch{iPatch}, 3).element_size);
                l4 = sum(msh_eval_boundary_side(obj.Meshes.msh_patch{iPatch}, 4).element_size);
                obj.Geometry(iPatch).Lengths = {l1; l2; l3; l4};
            end
            modelArea = sum(vertcat(obj.Geometry.Area));
            for iPatch = 1:obj.NumberPatches
                L1 = obj.Geometry(iPatch).Lengths{1} + obj.Geometry(iPatch).Lengths{2};
                L2 = obj.Geometry(iPatch).Lengths{3} + obj.Geometry(iPatch).Lengths{4};
                num_quivers_i = ceil(obj.Geometry(iPatch).Area/modelArea*nQuivers);
                % Number quivers in x an y according to lengths of sides
                q1 = ceil((num_quivers_i* L2/L1)^0.5);
                q2 = ceil((num_quivers_i* L1/L2)^0.5);
                % Points with also half distance to boundary
                paramPlotPoints = {linspace(0, 1-1/q1, q1)+0.5/q1, linspace(0, 1-1/q2, q2)+0.5/q2};
                % Points also directly on boundary
                %paramPlotPoints = {linspace(0, 1, q1), linspace(0, 1, q2)};
                obj.Geometry(iPatch).PlotPointsQuiverParam = paramPlotPoints;
                obj.Geometry(iPatch).PlotPointsQuiver = obj.getGeoPlotPoints(obj.Geometry(iPatch).nurbs, obj.Geometry(iPatch).PlotPointsQuiverParam, obj.RotationAngle);
            end
        end

        function calculatePlotValuesQuiver(obj)
            progressbar([obj.Name ': Calculating function space at quiver plot points: '])
            for iPatch = 1:obj.NumberPatches
                obj.Geometry(iPatch).PlotValuesQuiver = shp_eval1(obj.Spaces.sp_patch{iPatch}, obj.Geometry(iPatch), obj.Geometry(iPatch).PlotPointsQuiverParam, {'value', 'gradient'});
                progressbar(iPatch/obj.NumberPatches);
            end
            progressbar('finished');
        end

        function updatePlotPoints(obj, rot)
            obj.RotationAngle = rot;
            updatePlotPoints@IGA_Region(obj, rot)
            for iPatch = 1:obj.NumberPatches
                obj.Geometry(iPatch).PlotPointsQuiver = obj.getGeoPlotPoints(obj.Geometry(iPatch).nurbs, obj.Geometry(iPatch).PlotPointsQuiverParam, obj.RotationAngle);
            end
        end

        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.NumberDOF, obj.NumberDOF);
            if nargout == 2
                Mlin = sparse(obj.NumberDOF, obj.NumberDOF);
                for iMat = 1:numel(obj.Materials)
                    sigma = obj.Materials(iMat).Material.getSigma();
                    if isprop(obj.Materials(iMat).Material, "LamThickness")
                        d = obj.Materials(iMat).Material.LamThickness;
                        coeff = 1/12*sigma*d^2;
                        Mlin = Mlin + 1/obj.Length*op_gradu_gradv_mp_eval(obj.Spaces, obj.SpacesEval, obj.Spaces, obj.SpacesEval, obj.Meshes, obj.MeshesEval, coeff, obj.Materials(iMat).Patches);

                    else
                        Mlin = Mlin + 1/obj.Length*op_u_v_mp_eval(obj.Spaces, obj.SpacesEval, obj.Spaces, obj.SpacesEval, obj.Meshes, obj.MeshesEval, sigma, obj.Materials(iMat).Patches);
                    end
                end
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(obj.NumberDOF, obj.NumberDOF);
            for iMat = 1:numel(obj.Materials)
                if obj.Materials(iMat).Material.IsLinearMAG == false
                    Knonlin = Knonlin - 1/obj.Length*op_gradu_nu_gradv_mp_eval(obj.Spaces, obj.SpacesEval, obj.Spaces, obj.SpacesEval, ...
                        obj.Meshes, obj.MeshesEval, y/obj.Length, obj.Materials(iMat).Material , obj.Materials(iMat).Patches);
                end
            end
            if nargout == 2
                Klin = sparse(obj.NumberDOF, obj.NumberDOF);
                for iMat = 1:numel(obj.Materials)
                    if obj.Materials(iMat).Material.IsLinearMAG == true
                        nu = obj.Materials(iMat).Material.getNuLinear();
                        Klin = Klin - 1/obj.Length*op_gradu_gradv_mp_eval(obj.Spaces, obj.SpacesEval, obj.Spaces, obj.SpacesEval, obj.Meshes, obj.MeshesEval, nu, obj.Materials(iMat).Patches);
                    end
                end
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(obj.NumberDOF, 1);
            for iEXC = 1:numel(obj.Excitations)
                if ~obj.Excitations(iEXC).Excitation.IsConstant
                    Fnonlin = Fnonlin + obj.Excitations(iEXC).Excitation.getEXCvalues(t);
                end
            end
            if nargout == 2
                Flin = sparse(obj.NumberDOF, 1);
                for iEXC = 1:numel(obj.Excitations)
                    if obj.Excitations(iEXC).Excitation.IsConstant
                        Flin = Flin + obj.Excitations(iEXC).Excitation.getEXCvalues(t);
                    end
                end
            end
        end

        function J = JacobianMatrix(obj, t, y)
            J = sparse(obj.NumberDOF, obj.NumberDOF);
            for iMat = 1:numel(obj.Materials)
                if obj.Materials(iMat).Material.IsLinearMAG == false
                    J = J - 1/obj.Length*op_dK_du_times_u_mp_eval(obj.Spaces, obj.SpacesEval, obj.Spaces, obj.SpacesEval,...
                        obj.Meshes, obj.MeshesEval, y/obj.Length, obj.Materials(iMat).Material, obj.Materials(iMat).Patches);
                end
            end
        end

        function postprocess(obj, ts, ys)
            obj.Solution.Time = ts;
            obj.Solution.Values = ys / obj.Length;
            if isscalar(ts)
                obj.MagneticPotential = ys / obj.Length;
            end
            if ~isempty(obj.EvaluationShapeValues)
                obj.EvaluationValues = squeeze(sum((sum(double(obj.EvaluationShapeValues{2}).*permute(ys/obj.Length, [4,3,1,2]), 3)).^2, 2).^0.5);
            end
        end

        %%%%%%%%%%%% PLOTTING

        function checkPlotValues(obj)
            if obj.PlotValuesNeedUpdate
                obj.calculatePlotValues();
                obj.calculatePlotValuesQuiver();
            end
            obj.PlotValuesNeedUpdate = false;
        end

        function plotMagneticFluxDensity(obj, t, cols)
            if ~exist("t", "var")
                t = 0;
            end
            if size(obj.Solution.Values, 2) > 1
                obj.MagneticPotential = interp1(obj.Solution.Time, obj.Solution.Values', t)';
            else
                obj.MagneticPotential = obj.Solution.Values;
            end
            warning("off", "MATLAB:contour:ConstantData")
            if ~exist('cols', 'var')
                cols = [0, 2];
            end
            obj.checkPlotValues()

            % if figure already exists, update current one
            if ~ishandle(obj.FigureBRes)
                obj.FigureBRes = figure('Name', 'B Resulting Plot', 'NumberTitle', 'off', 'Position', [0 0 800 500]);
                hold on
                c = colorbar();
                c.TickLabelInterpreter = 'latex';
                c.Label.Interpreter = 'latex';
                c.Label.String = 'Magnetic flux density (T)';
                c.FontSize = 18;
                c.Label.FontSize = 20;
                % set(c, "Visible", "off");

                view([0, 90]);
                axis equal;
            else
                figure(obj.FigureBRes); % focus existing figure
            end
            % Plot contour
            for iPatch = 1:obj.NumberPatches
                a_i = obj.MagneticPotential(obj.Spaces.gnum{iPatch});
                B_xy = squeeze(sum(sum(obj.Geometry(iPatch).PlotValues{2}.* reshape(a_i, 1, 1, 1, []), 4).^2, 1).^0.5);
                if ishandle(obj.FigureBResBMag(iPatch))
                    set(obj.FigureBResBMag(iPatch), 'XData', obj.Geometry(iPatch).PlotPoints{1}, 'YData', obj.Geometry(iPatch).PlotPoints{2}, 'ZData', B_xy, 'LevelList', linspace(cols(1), cols(2), 20));
                else
                    [~, obj.FigureBResBMag(iPatch)] = contourf(obj.Geometry(iPatch).PlotPoints{1}, obj.Geometry(iPatch).PlotPoints{2}, B_xy, linspace(cols(1), cols(2), 20));
                    set(obj.FigureBResBMag(iPatch), 'LineColor', 'none');
                end
            end
            % Plot interfaces
            for iInt = 1:numel(obj.Interfaces)
                if ishandle(obj.FigureBResInterfaces(iInt))
                    set(obj.FigureBResInterfaces(iInt), 'XData', obj.Interfaces(iInt).PlotPoints{1}, 'YData', obj.Interfaces(iInt).PlotPoints{2});
                else
                    obj.FigureBResInterfaces(iInt) = plot(obj.Interfaces(iInt).PlotPoints{1}, obj.Interfaces(iInt).PlotPoints{2}, 'Color', obj.Interfaces(iInt).PlotColor);
                end
            end
            % Plot boundaries
            for iBnd = 1:numel(obj.Boundaries)
                if ishandle(obj.FigureBResBoundaries(iBnd))
                    set(obj.FigureBResBoundaries(iBnd), 'XData', obj.Boundaries(iBnd).PlotPoints{1}, 'YData', obj.Boundaries(iBnd).PlotPoints{2});
                else
                    obj.FigureBResBoundaries(iBnd) = plot(obj.Boundaries(iBnd).PlotPoints{1}, obj.Boundaries(iBnd).PlotPoints{2}, 'Color', obj.Boundaries(iBnd).PlotColor);
                end
            end
            % Plot quivers
            X = cell(obj.NumberPatches, 1);
            Y = cell(obj.NumberPatches, 1);
            BX = cell(obj.NumberPatches, 1);
            BY = cell(obj.NumberPatches, 1);
            for iPatch = 1:obj.NumberPatches
                a_i = obj.MagneticPotential(obj.Spaces.gnum{iPatch});
                Bx =  sum(obj.Geometry(iPatch).PlotValuesQuiver{2}(2, :, :, :).* reshape(a_i, 1, 1, 1, []), 4);
                By = -sum(obj.Geometry(iPatch).PlotValuesQuiver{2}(1, :, :, :).* reshape(a_i, 1, 1, 1, []), 4);
                X{iPatch} = reshape(obj.Geometry(iPatch).PlotPointsQuiver{1, 1}, [], 1);
                Y{iPatch} = reshape(obj.Geometry(iPatch).PlotPointsQuiver{2, 1}, [], 1);
                BX{iPatch} = reshape(Bx, [], 1);
                BY{iPatch} = reshape(By, [], 1);
            end
            X = cell2mat(X)';
            Y = cell2mat(Y)';
            BX = cell2mat(BX)';
            BY = cell2mat(BY)';
            % Rotation
            Bxtemp = BX*cos(obj.RotationAngle) - BY*sin(obj.RotationAngle);
            Bytemp = BY*cos(obj.RotationAngle) + BX*sin(obj.RotationAngle);
            BX = Bxtemp;
            BY = Bytemp;

            if ishandle(obj.FigureBResQuiver)
                set(obj.FigureBResQuiver, 'XData', X, 'YData', Y, 'UData', BX, 'VData', BY);
            else
                obj.FigureBResQuiver = quiver(X, Y, BX, BY, 'color', 'k');
            end
            xlim(obj.PlotLimitsX);
            ylim(obj.PlotLimitsY);
            clim(cols);
            drawnow();
        end

        function plotMagneticPotentialLines(obj, t, nlines, limits)
            if ~exist("t", "var")
                t = 0;
            end
            if size(obj.Solution.Values, 2) > 1
                obj.MagneticPotential = interp1(obj.Solution.Time, obj.Solution.Values', t)';
            else
                obj.MagneticPotential = obj.Solution.Values;
            end
            if ~exist('nlines', 'var')
                nlines = 20;
            elseif isempty(nlines)
                nlines = 20;
            end
            if ~exist('limits', 'var')
                limits = [min(obj.MagneticPotential), max(obj.MagneticPotential)];
            end
            obj.checkPlotValues();

            if ~ishandle(obj.FigurePotLines)
                % create new plot for Geometry
                obj.FigurePotLines = figure('Name', 'Potential Line Plot', 'NumberTitle', 'off', 'Position', [0 0 800 500]);
                hold on
                axis equal
                view([0, 90]);
            else
                figure(obj.FigurePotLines);
            end
            % Potential Lines
            for iPatch = 1:obj.NumberPatches
                a_i = obj.MagneticPotential(obj.Spaces.gnum{iPatch});
                Az = obj.Geometry(iPatch).PlotValues{1}.* reshape(a_i, 1, 1, []);
                Az = sum(Az, 3);
                if ishandle(obj.FigurePotLinesGeometry(iPatch))
                    set(obj.FigurePotLinesGeometry(iPatch), 'XData', obj.Geometry(iPatch).PlotPoints{1}, 'YData', obj.Geometry(iPatch).PlotPoints{2});
                else
                    iMat = find(cellfun(@(x) any(x==iPatch), {obj.Materials.Patches}));
                    obj.FigurePotLinesGeometry(iPatch)=surf(obj.Geometry(iPatch).PlotPoints{1}, obj.Geometry(iPatch).PlotPoints{2}, zeros(size(obj.Geometry(iPatch).PlotPoints{1})), ...
                        "EdgeColor", "none", "FaceColor", obj.Materials(iMat).Material.getPlotColor(), "FaceAlpha", obj.Materials(iMat).Material.getPlotAlpha());
                end
                if ishandle(obj.FigurePotLinesPotentialLines(iPatch))
                    set(obj.FigurePotLinesPotentialLines(iPatch), 'XData', obj.Geometry(iPatch).PlotPoints{1}, 'YData', obj.Geometry(iPatch).PlotPoints{2}, 'ZData', Az, 'LevelList', linspace(limits(1), limits(2), nlines));
                else
                    [~, obj.FigurePotLinesPotentialLines(iPatch)] = contour(obj.Geometry(iPatch).PlotPoints{1}, obj.Geometry(iPatch).PlotPoints{2}, Az, linspace(limits(1), limits(2), nlines), LineColor="black");%, LineColor="black");
                end
            end
            % Interfaces
            for iInt = 1:numel(obj.Interfaces)
                if ishandle(obj.FigurePotLinesInterfaces(iInt))
                    set(obj.FigurePotLinesInterfaces(iInt), 'XData', obj.Interfaces(iInt).PlotPoints{1}, 'YData', obj.Interfaces(iInt).PlotPoints{2});
                else
                    obj.FigurePotLinesInterfaces(iInt) = plot(obj.Interfaces(iInt).PlotPoints{1}, obj.Interfaces(iInt).PlotPoints{2}, "Color", obj.Interfaces(iInt).PlotColor);
                end
            end
            % Boundaries
            for iBnd = 1:numel(obj.Boundaries)
                if ishandle(obj.FigurePotLinesBoundaries(iBnd))
                    set(obj.FigurePotLinesBoundaries(iBnd), 'XData', obj.Boundaries(iBnd).PlotPoints{1}, 'YData', obj.Boundaries(iBnd).PlotPoints{2});
                else
                    obj.FigurePotLinesBoundaries(iBnd) = plot(obj.Boundaries(iBnd).PlotPoints{1}, obj.Boundaries(iBnd).PlotPoints{2}, "Color", obj.Boundaries(iBnd).PlotColor);
                end
            end
            drawnow();
            %shading interp
            colormap jet
            xlim(obj.PlotLimitsX);
            ylim(obj.PlotLimitsY);
        end

        function plotMagneticVectorPotential(obj)
            obj.checkPlotValues();
            % figure
            hold on
            for iPatch = 1:obj.NumberPatches
                a_i = obj.MagneticPotential(obj.Spaces.gnum{iPatch});
                Az = sum(obj.Geometry(iPatch).PlotValues{1}.*reshape(a_i, 1, 1, []), 3);
                surf(obj.Geometry(iPatch).PlotPoints{1}, obj.Geometry(iPatch).PlotPoints{2}, Az);
            end
            % boundaries
            % for iBnd = 1:numel(obj.Boundaries)
            %     obj.FigureBResBoundaries(iBnd) = plot(obj.Boundaries(iBnd).PlotPoints{1}, obj.Boundaries(iBnd).PlotPoints{2}, "Color", obj.Boundaries(iBnd).PlotColor);
            % end
            % % interfaces
            % for iInt = 1:numel(obj.Interfaces)
            %     obj.FigureBResInterfaces(iInt) = plot(obj.Interfaces(iInt).PlotPoints{1}, obj.Interfaces(iInt).PlotPoints{2}, "Color", obj.Interfaces(iInt).PlotColor);
            % end

            view([0, 90]);
            % axis equal;
            xlim(obj.PlotLimitsX);
            ylim(obj.PlotLimitsY);
        end

        function initializeFigures(obj)
            initializeFigures@IGA_Region(obj);
            obj.FigureBRes = gobjects(1);
            obj.FigureBResInterfaces = gobjects(numel(obj.Interfaces), 1);
            obj.FigureBResBoundaries = gobjects(numel(obj.Boundaries), 1);
            obj.FigureBResBMag = gobjects(obj.NumberPatches, 1);
            obj.FigureBResQuiver = gobjects(obj.NumberPatches, 1);

            obj.FigurePotLines = gobjects(1);
            obj.FigurePotLinesInterfaces = gobjects(numel(obj.Interfaces), 1);
            obj.FigurePotLinesBoundaries = gobjects(numel(obj.Boundaries), 1);
            obj.FigurePotLinesGeometry = gobjects(obj.NumberPatches, 1);
            obj.FigurePotLinesPotentialLines = gobjects(obj.NumberPatches, 1);

            obj.FigureRemanence = gobjects(1);
        end

        function exportParaviewAz(obj, filename, valueIndex, valueList)
            if ~isfolder("paraview")
                mkdir("paraview")
            end
            if ~exist("filename", "var")
                filename = 'paraview/Az';
            end
            if nargin <= 2
                valueIndex = 1;
                valueList = 1;
            end
            obj.checkPlotValues();

            str1 = cat (2, '<?xml version="1.0"?> \n', '<VTKFile type="Collection" version="0.1"> \n', '<Collection> \n');
            str2 = cat (2, '<DataSet timestep="%d" group="" part="%d" file="%s.vts"/> \n');
            str3 = cat (2, '</Collection>\n', '</VTKFile> \n');

            if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.pvd'))
                pvd_filename = cat (2, filename, '.pvd');
            else
                pvd_filename = filename;
                filename = filename (1:end-4);
            end

            if valueIndex == 1
                fid = fopen (pvd_filename, 'w');
                fprintf (fid, str1);
            else
                fid = fopen (pvd_filename, 'a');
            end

            if (fid < 0)
                error ('exportSolutionParaView: could not open file %s', pvd_filename);
            end

            ind = union (find (filename == '/', 1, 'last'), find (filename == '\', 1, 'last')) + 1;
            if (isempty (ind)); ind = 1; end

            for iPatch = 1:obj.NumberPatches
                filename_patch_without_path = cat (2, filename(ind:end), '_', num2str(valueIndex-1), '_', num2str (iPatch));
                filename_patch = cat (2, filename, '_', num2str(valueIndex-1), '_', num2str (iPatch));
                vts_pts = [];

                vts_pts(1, :, :) = obj.Geometry(iPatch).PlotPoints{1};
                vts_pts(2, :, :) = obj.Geometry(iPatch).PlotPoints{2};
                Az_patch = obj.MagneticPotential(obj.Spaces.gnum{iPatch});
                Az_points = sum(obj.Geometry(iPatch).PlotValues{1}.* reshape(Az_patch, 1, 1, []), 3);

                msh_to_vtk (vts_pts, Az_points, filename_patch, 'A_z');

                fprintf (fid, str2, valueList(valueIndex), iPatch, filename_patch_without_path);
            end

            if valueIndex == numel(valueList)
                fprintf (fid, str3);
            end
            fclose (fid);
        end

        function exportTransientParaview(obj, filename, times)
            progressbar('Exporting to Paraview: ')
            for ti = 1:numel(times)
                obj.MagneticPotential = interp1(obj.Solution.Time, obj.Solution.Values', times(ti))';
                obj.exportParaviewAz(['paraview/' filename '_Az'], ti, times);
                obj.exportParaviewB(['paraview/' filename '_B'], ti, times);
                progressbar(ti/numel(times));
            end
            progressbar('finished');
        end

        %% Export Magnetic B Field to ParaView
        function exportParaviewB(obj, filename, currentIndex, dataList)
            if ~isfolder("paraview")
                mkdir("paraview")
            end
            if ~exist("filename", "var")
                filename = 'paraview/B';
            end
            if nargin <= 2
                currentIndex = 1;
                dataList = 1;
            end
            obj.checkPlotValues();

            str1 = cat (2, '<?xml version="1.0"?> \n', '<VTKFile type="Collection" version="0.1"> \n', '<Collection> \n');
            str2 = cat (2, '<DataSet timestep="%d" group="" part="%d" file="%s.vts"/> \n');
            str3 = cat (2, '</Collection>\n', '</VTKFile> \n');

            if (length (filename) < 4 || ~strcmp (filename(end-3:end), '.pvd'))
                pvd_filename = cat (2, filename, '.pvd');
            else
                pvd_filename = filename;
                filename = filename (1:end-4);
            end

            if currentIndex == 1
                fid = fopen (pvd_filename, 'w');
                fprintf (fid, str1);
            else
                fid = fopen (pvd_filename, 'a');
            end

            if (fid < 0)
                error ('exportSolutionParaView: could not open file %s', pvd_filename);
            end

            ind = union (find (filename == '/', 1, 'last'), find (filename == '\', 1, 'last')) + 1;
            if (isempty (ind)); ind = 1; end

            for iPatch = 1:obj.NumberPatches
                filename_patch_without_path = cat (2, filename(ind:end), '_', num2str(currentIndex-1), '_', num2str (iPatch));
                filename_patch = cat (2, filename, '_', num2str(currentIndex-1), '_', num2str (iPatch));
                vts_pts = [];

                vts_pts(1, :, :) = obj.Geometry(iPatch).PlotPointsQuiver{1,1};
                vts_pts(2, :, :) = obj.Geometry(iPatch).PlotPointsQuiver{2,1};
                Az_patch = obj.MagneticPotential(obj.Spaces.gnum{iPatch});

                Bx = sum(obj.Geometry(iPatch).PlotValuesQuiver{2}(2, :, :, :).* reshape(Az_patch, 1, 1, 1, []), 4);
                By = -sum(obj.Geometry(iPatch).PlotValuesQuiver{2}(1, :, :, :).* reshape(Az_patch, 1, 1, 1, []), 4);
                angle = 0;
                Bxtemp = Bx*cos(angle) - By*sin(angle);
                Bytemp = By*cos(angle) + Bx*sin(angle);
                B = [Bxtemp; Bytemp];

                msh_to_vtk_vec (vts_pts, B, filename_patch, 'B');
                fprintf (fid, str2, dataList(currentIndex), iPatch, filename_patch_without_path);
            end

            if currentIndex == numel(dataList)
                fprintf (fid, str3);
            end
            fclose (fid);
        end
    end
end