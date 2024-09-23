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
Rotor = IGA_RegionMagnetic_OPT();
Rotor.importSurface(srfRotor);
Rotor.setProperties([4, 4], [2,2], 0.035);
Rotor.calculatePlotPoints([5, 5]);

Stator = IGA_RegionMagnetic_OPT();
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
MatMagnetLeft =  MAT_NdFeB_Br10(MA2);
% Redefine Cost/kg, such that MagnetArea = Cost
Iron.Cost=0;
MatMagnetRight.Cost = 1/(MatMagnetRight.Rho*0.035);
MatMagnetLeft.Cost = 1/(MatMagnetLeft.Rho*0.035);

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
PhaseU = IGA_EXCcurrentSine_OPT(Stator, phaseDefU, I0, pp, 0,      phase0, NumberWindings);
PhaseV = IGA_EXCcurrentSine_OPT(Stator, phaseDefV, I0, pp,-2/3*pi, phase0, NumberWindings);
PhaseW = IGA_EXCcurrentSine_OPT(Stator, phaseDefW, I0, pp, 2/3*pi, phase0, NumberWindings);

Stator.resetExcitations()
Stator.addExcitation(PhaseU);
Stator.addExcitation(PhaseV);
Stator.addExcitation(PhaseW);

%% Coupling
Coupling = IGA_HarmonicMortaring_OPT(Rotor, [2:5, 10, 12:15], Stator, [2, 66, 7:2:63, 71:2:79]);

nharm = 26;
mult = opts.POLES;
% Sinus and Cosinus values: 2, 6, 10, 14, ...
SinValues = mult/2:mult:mult*nharm;
CosValues = mult/2:mult:mult*nharm;
Coupling.setHarmonics(SinValues, CosValues);

%%
Rotor.resetOptimizationResponses();
Rotor.resetOptimizationParameters();
Rotor.addOptimizationResponse("Cost");
Rotor.addOptimizationResponse("SurfaceSmoothness");
Rotor.addGeometryConstraintFunction("RotorConstraints", @(p) GEO_IPMrotor_V_1_constr(p));

optiPoints =    {367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 354, 394, 395, 396, 397, 392, 403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 331, 334};
optiPointsSym = {246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 233, 273, 274, 275, 276, 271, 284, 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, 282, 333};

% optiPoints =    {[367, 368, 369], 370, 371, 372, 373, 374, 375, [376, 354, 394], 395, 396, [397, 392, 403], 404, 405, 406, 407, 408, 409, 410, 411, 412, [413, 331, 334]};
% optiPointsSym = {[246, 247, 248], 249, 250, 251, 252, 253, 254, [255, 233, 273], 274, 275, [276, 271, 284], 285, 286, 287, 288, 289, 290, 291, 292, 293, [294, 282, 333]};

for iPoint = 1:numel(optiPoints)
    Rotor.addOptimizationControlPointOffset("CP"+num2str(iPoint), [optiPoints{iPoint}, optiPointsSym{iPoint}],"r", 0, -1.5e-3, 0)
end

Rotor.addOptimizationControlPointOffset("CP"+num2str(numel(optiPoints)+1), [344, 223], "r", 0, -1.0e-3, 0);
Rotor.addOptimizationControlPointOffset("CP"+num2str(numel(optiPoints)+2), [366, 245], "r", 0, -1.0e-3, 0);
Rotor.addOptimizationControlPointOffset("CP"+num2str(numel(optiPoints)+3), [348, 227], "r", 0, -1.0e-4, 0);
%%
Stator.resetOptimizationResponses();
Stator.resetOptimizationParameters();

Stator.addOptimizationParameter("Phase0", 0, -deg2rad(20), deg2rad(20));

Rotor.GeometryFunction = @(p) GEO_IPMrotor_V_1(p);

Coupling.resetOptimizationResponses();
Coupling.addOptimizationResponse("Torque");
Coupling.addOptimizationResponse("TorqueSTD");
Coupling.addOptimizationResponse("TorqueMEAN");

%%
Motor = IGA_MotorClass_OPT(Rotor, Stator);
Motor.setCoupling(Coupling);

Motor.OptimizationAngles = deg2rad([0:1:29]);

%%
Motor.resetOptimizationConstraints();
Motor.resetOptimizationObjectives();

Motor.addOptimizationObjective("Cost", 1e4);
Motor.addOptimizationObjective("TorqueSTD", 100);
Motor.addOptimizationObjective("SurfaceSmoothness", 1e3);

Motor.addOptimizationConstraint("TorqueMEAN", ">", 2/4, 100);

Motor.addOptimizationConstraint("RotorConstraints1", "<", 0)
Motor.addOptimizationConstraint("RotorConstraints2", "<", 0)
Motor.addOptimizationConstraint("RotorConstraints3", "<", 0)
Motor.addOptimizationConstraint("RotorConstraints4", "<", 0)
Motor.addOptimizationConstraint("RotorConstraints5", "<", 0)
Motor.addOptimizationConstraint("RotorConstraints6", "<", 0)
Motor.addOptimizationConstraint("RotorConstraints7", "<", 0)
Motor.addOptimizationConstraint("RotorConstraints8", "<", 0)
Motor.addOptimizationConstraint("RotorConstraints9", "<", 0)


%%
xopt = Motor.optimize();
