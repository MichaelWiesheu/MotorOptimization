clear all
close all
clc
restoredefaultpath

addpath(genpath("../../"))

%%
folder = '../opt_param/plots/final';
OptimizationHistory = load([folder '/OptimizationHistory']).OptimizationHistory;
opts.MA = OptimizationHistory(end).Parameters(8).Value;

%%
opts.draw_geometry = false;
opts.POLES = 4;

folder = '../opt_param_then_shape/plots/final';
srfName = '/RotorOpt220.txt';
OptimizationHistory = load([folder '/OptimizationHistory']).OptimizationHistory;

opts.Phase0 = OptimizationHistory(end).Parameters(33).Value;

[srfRotor, patchesRotor, opts] = GEO_IPMrotor_V_1(opts);
srfRotor = [folder srfName];

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
phase0 = opts.Phase0;
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
xlim([0, 0.055]);
ylim([0, 0.055]);
axis off
ax = gca;
% Requires R2020a or later
exportgraphics(ax,['plots/paramthenshapeOpt.pdf']) 
exportgraphics(ax,['plots/paramthenshapeOpt.png'], 'Resolution', 300) 

%%
y0 = zeros(Motor.NumberDOF, 1);
angles = 0:0.25:30;
clear T T1 T2
for iAngle = 1:numel(angles)
    angle = angles(iAngle);
    Motor.setRotationAngle(angle, "deg");
    y0 = Motor.solveStatic(deg2rad(angle), y0);
    T(iAngle) = Coupling.Solution.Torque;
    T1(iAngle) = Coupling.calcTorqueBrBtRt();
    T2(iAngle) = Coupling.calcTorqueBrBtSt();
end
save("data/Torque_param_then_shape.mat", "angles", "T", "T1", "T2")

