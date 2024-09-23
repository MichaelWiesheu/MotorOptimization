classdef IGA_EXCcurrent < EXC
    properties (Access = public)
        RHS;
        Current;
        Phases;
        NumberWindings;
        Patches;
        Color = [0,0,0];
    end

    methods
        function obj = IGA_EXCcurrent(parentRegion, phasesStruct, currentFunction, nWindings)
            obj@EXC(parentRegion);
            obj.Current = currentFunction;
            if isa(currentFunction, "function_handle")
                obj.IsConstant = false;
            else
                obj.IsConstant = true;
            end
            obj.NumberWindings = nWindings;
            for iPhase = 1:numel(phasesStruct)
                obj.Phases(iPhase).Patches = phasesStruct(iPhase).Patches;
                obj.Phases(iPhase).Slot = phasesStruct(iPhase).Slot;
            end
            obj.Patches = horzcat(phasesStruct.Patches);
            obj.updateExcitation()
        end

        function updateExcitation(obj)
            obj.RHS = zeros(obj.ParentRegion.NumberDOF, 1);
            for iPhase = 1:numel(obj.Phases)
                patches = obj.Phases(iPhase).Patches;
                slot = obj.Phases(iPhase).Slot;
                obj.Phases(iPhase).Area = sum(vertcat(obj.ParentRegion.Geometry(patches).Area));
                obj.Phases(iPhase).RHS = op_f_v_mp_eval (obj.ParentRegion.Spaces, obj.ParentRegion.SpacesEval, obj.ParentRegion.MeshesEval, patches);
                obj.RHS = obj.RHS + slot*obj.NumberWindings/obj.Phases(iPhase).Area* obj.Phases(iPhase).RHS;
            end
        end

        function f = getEXCvalues(obj, t)
            if obj.IsConstant
                f = obj.Current*obj.RHS;
            else
                f = obj.Current(t)*obj.RHS;
            end
        end

        function plotWindingCurrent(obj, radialPlot)
            if ~exist("radialPlot", "var")
                radialPlot = true;
            end
            
            for iPhase = 1:numel(obj.Phases)
                AreaMoment = 0;
                for iPatch = 1:numel(obj.Phases(iPhase).Patches)
                    patch = obj.Phases(iPhase).Patches(iPatch);
                    AreaMoment = AreaMoment + obj.ParentRegion.Geometry(patch).Area*nrbeval(obj.ParentRegion.Geometry(patch).nurbs, {0.5, 0.5});
                end
                Xmean = AreaMoment/obj.Phases(iPhase).Area;
                if radialPlot
                    angle = atan2(Xmean(2), Xmean(1));
                else
                    angle = 0;
                end
                if obj.Phases(iPhase).Slot == 1
                    plotInCurrent(Xmean(1), Xmean(2), obj.Phases(iPhase).Area^0.5/5);
                else
                    plotOutCurrent(Xmean(1), Xmean(2), obj.Phases(iPhase).Area^0.5/5, angle);
                end
            end

            function plotInCurrent(x, y, r)
                p = nsidedpoly(1000, 'Center', [x y], 'Radius', r);
                p1 = nsidedpoly(1000, 'Center', [x y], 'Radius', r/5);
                plot(p, 'FaceColor', 'white', 'EdgeColor', obj.Color, 'FaceAlpha', 1, 'LineWidth', 3)
                plot(p1, 'FaceColor', obj.Color, 'EdgeColor', obj.Color, 'FaceAlpha', 1)
            end

            function plotOutCurrent(x, y, r, rot)
                p = nsidedpoly(1000, 'Center', [x y], 'Radius', r);
                pts = [-r, -r, r, r; -r, r -r, r]/(2^0.5);
                R = [cos(rot), -sin(rot); sin(rot), cos(rot)];
                pts = R*pts;
                plot(p, 'FaceColor', 'white', 'EdgeColor', obj.Color, 'FaceAlpha', 1, 'LineWidth', 3)
                plot(x + [pts(1, 1), pts(1, 4) ], y + [pts(2, 1), pts(2, 4) ], 'LineWidth', 3, 'Color', obj.Color)
                plot(x + [pts(1, 2), pts(1, 3) ], y + [pts(2, 2), pts(2, 3) ], 'LineWidth', 3, 'Color', obj.Color)
            end
        end
    end
end