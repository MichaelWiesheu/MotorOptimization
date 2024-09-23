classdef IGA_MotorClass < DAE_System
    
    properties
        Rotor, Stator;
        Coupling;
        FigureGeometry
        PlotLimitsX=[-Inf,Inf], PlotLimitsY=[-Inf,Inf];
    end
    
    methods
        function obj = IGA_MotorClass(rotor, stator)
            obj.Rotor = rotor;
            obj.Stator = stator;
            obj.addElement(rotor);
            obj.addElement(stator);
        end

        function setCoupling(obj, coupling)
            obj.Coupling = coupling;
            obj.addElement(coupling)
        end

        function setRotationAngle(obj, angle, unit)
            if ~exist("unit", "var")
                unit = "rad";
            end
            if unit == "deg"
                angle = deg2rad(angle);
            elseif unit == "rad"
                % pass
            else
                warning("Unit not known");
            end
            obj.Coupling.setRotationAngle(angle);
        end

        function plotGeometry(obj)
            obj.Rotor.plotGeometry()
            obj.Stator.FigureGeometry = obj.Rotor.FigureGeometry;
            obj.Stator.plotGeometry()
        end

        function plotMagneticFluxDensity(obj, t)
            if ~exist("t", "var")
                t = 0;
            end
            obj.Rotor.plotMagneticFluxDensity(t)
            obj.Stator.FigureBRes = obj.Rotor.FigureBRes;
            obj.Stator.plotMagneticFluxDensity(t)
        end

        function plotMagneticPotentialLines(obj, t, nlines)
            if ~exist('t', 'var')
                t = 0;
            end
            if ~exist('nlines', 'var')
                nlines = 20;
            end
            limits = [min([obj.Rotor.MagneticPotential; obj.Stator.MagneticPotential]), ...
                    max([obj.Rotor.MagneticPotential; obj.Stator.MagneticPotential])];

            obj.Rotor.plotMagneticPotentialLines(t, nlines, limits)
            obj.Stator.FigurePotLines = obj.Rotor.FigurePotLines;
            obj.Stator.plotMagneticPotentialLines(t, nlines, limits)
        end
    end
end
