
clear all
close all
clc
restoredefaultpath

addpath(genpath("../../"))

%%
opts.draw_geometry = false;
opts.MA = 150;
opts.POLES = 4;

[srfRotor, patchesRotor, opts] = GEO_IPMrotor_V_1(opts);
[srfStator, patchesStator] = GEO_DISstator_2(opts);

%%
Rotor = IGA_RegionMagnetic();
Rotor.importSurface(srfRotor);
Rotor.setProperties([4, 4], [2,2], 0.035);
Rotor.calculatePlotPoints([5, 5]);

Stator = IGA_RegionMagnetic();
Stator.importSurface(srfStator);
Stator.setProperties([3, 3], [2,2], 0.035);
Stator.calculatePlotPoints([5,5]);

%% Materials
Air = MAT_Air();
Iron = MAT_M27();
Copper = MAT_Copper();
MA1 = pi/2-(deg2rad(opts.MA/2) - pi/4);
MA2 = MA1 -pi + deg2rad(opts.MA);

MatMagnetRight = MAT_NdFeB_Br10(MA1);
MatMagnetLeft = MAT_NdFeB_Br10(MA2);

Rotor.resetMaterials();
Rotor.setMaterial(Air, patchesRotor.Air);
Rotor.setMaterial(Iron, patchesRotor.Iron);
Rotor.setMaterial(MatMagnetRight, patchesRotor.Magnets(1));
Rotor.setMaterial(MatMagnetLeft, patchesRotor.Magnets(2));

Stator.resetMaterials();
Stator.setMaterial(Air, patchesStator.Air);
Stator.setMaterial(Iron, patchesStator.Iron);
Stator.setMaterial(Copper, patchesStator.Windings);

%%
Rotor.resetBoundaryConditions()
Rotor.addBoundaryCondition(IGA_BCdirichlet(Rotor, [8, 9, 18]));
Rotor.addBoundaryCondition(IGA_BCantiPeriodic(Rotor, [16, 17, 11], [6, 7, 1]));

%%
Stator.resetBoundaryConditions()
Stator.addBoundaryCondition(IGA_BCdirichlet(Stator, [6:2:64, 70:2:80]));
Stator.addBoundaryCondition(IGA_BCantiPeriodic(Stator, [65, 67, 68, 69], [1, 3, 4, 5]));

%%
phaseDefU(1).Patches = [7, 11];         phaseDefU(1).Slot = 1;
phaseDefU(2).Patches = [19, 23];        phaseDefU(2).Slot = 1;
phaseDefU(3).Patches = [43, 47];        phaseDefU(3).Slot = 1;
phaseDefU(4).Patches = [31, 35];        phaseDefU(4).Slot = 1;

phaseDefW(1).Patches = [67, 71];        phaseDefW(1).Slot = -1;
phaseDefW(2).Patches = [55, 59];        phaseDefW(2).Slot = -1;
phaseDefW(3).Patches = [91, 95];        phaseDefW(3).Slot = -1;
phaseDefW(4).Patches = [79, 83];        phaseDefW(4).Slot = -1;

phaseDefV(1).Patches = [115, 119];      phaseDefV(1).Slot = 1;
phaseDefV(2).Patches = [103, 107];      phaseDefV(2).Slot = 1;
phaseDefV(3).Patches = [139, 143];      phaseDefV(3).Slot = 1;
phaseDefV(4).Patches = [127, 131];      phaseDefV(4).Slot = 1;


NumberWindings = 35;
phase0 = 0;
I0 = 3;
pp = 2;
% synchronous speed: omega*t = theta*pp
PhaseU = IGA_EXCcurrentSine(Stator, phaseDefU, I0, pp, 0,      phase0, NumberWindings);
PhaseV = IGA_EXCcurrentSine(Stator, phaseDefV, I0, pp,-2/3*pi, phase0, NumberWindings);
PhaseW = IGA_EXCcurrentSine(Stator, phaseDefW, I0, pp, 2/3*pi, phase0, NumberWindings);

Stator.resetExcitations()
Stator.addExcitation(PhaseU);
Stator.addExcitation(PhaseV);
Stator.addExcitation(PhaseW);

%% Coupling
Coupling = IGA_HarmonicMortaring(Rotor, [2:5, 10, 12:15], Stator, [2, 66, 7:2:63, 71:2:79]);

nharm = 26;
mult = opts.POLES;
% Sinus and Cosinus values: 2, 6, 10, 14, ...
SinValues = mult/2:mult:mult*nharm;
CosValues = mult/2:mult:mult*nharm;
Coupling.setHarmonics(SinValues, CosValues);


%%
Motor = IGA_MotorClass(Rotor, Stator);
Motor.setCoupling(Coupling);


%%
Motor.plotGeometry()
PhaseU.plotWindingCurrent();
PhaseV.plotWindingCurrent();
PhaseW.plotWindingCurrent();

xlim([0, Inf]);
ylim([0, Inf]);
axis off
% set(gca, "FontSize", 14)
% saveas(gcf, "Motor1.png")
%%
Motor.solveStatic();
Motor.plotMagneticFluxDensity();
axis off
set(gca, 'FontSize', 20)
ax = gca;
% Requires R2020a or later
exportgraphics(ax,'JMAGfield.pdf','Resolution',600) 
exportgraphics(ax,'JMAGfield.png') 
%%
y0 = zeros(Motor.NumberDOF, 1);
angles = 0:0.25:90;
clear T T1 T2
for iAngle = 1:numel(angles)
    angle = angles(iAngle);
    Motor.setRotationAngle(angle, "deg");
    y0 = Motor.solveStatic(deg2rad(angle), y0);
    T(iAngle) = Coupling.Solution.Torque;
    T1(iAngle) = Coupling.calcTorqueBrBtRt();
    T2(iAngle) = Coupling.calcTorqueBrBtSt();
end
save("Torque_init.mat", "angles", "T", "T1", "T2")
%%

set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(groot, 'defaultLegendInterpreter', 'latex');
% close f
f = figure('Position',[200,200,600,300]);
hold on

lw=2;
JMAG_torque05 = cell2mat(table2cell(readtable("JMAG_torque_init_05.txt")));
plot(JMAG_torque05(:, 1), JMAG_torque05(:, 2), LineWidth=lw, Color="black");

% scatter(angles, 4*T, "filled", "MarkerFaceColor", TUDa_getColor("7b"), SizeData=50)
plot(angles, 4*T, LineWidth=lw, Color=TUDa_getColor('9b'), LineStyle="-.");

ylim([1, 2.8])
xlim([0,90]);
legend(["JMAG", "Matlab"],Location="southeast")
grid on
xlabel("Rotation angle")
ylabel("Torque (Nm)")
set(gca, "FontSize", 15)

%%
jmagFun = @(val) interp1(JMAG_torque05(:, 1), JMAG_torque05(:, 2), val);

jmagTorque = jmagFun(angles);

rms(4*T-jmagTorque)
mean(abs(4*T-jmagTorque)./jmagTorque)
% mean 0.7% differences, explained by material inter/extrapolation or meshing.
 




