clc
clear all
close all
restoredefaultpath

addpath(genpath("../../"))

%%
JMAGoutput = readtable("data/JMAGOutput50");
%%
Case = JMAGoutput.Case(JMAGoutput.TargetTorque >= 2.3);
fopt = JMAGoutput.MinMagnetMass(JMAGoutput.TargetTorque >= 2.3);
%%
[MinVal, indx] = min(fopt);
Case(indx)
indxOpt = 2248;
%%
for i = 1:30
    eval("T("+string(i)+")=JMAGoutput.T"+string(i-1)+"(indxOpt);")
end
MT1 = JMAGoutput.CADParameters_MT1_Variables(indxOpt);
MW1 = JMAGoutput.CADParameters_MW1_Variables(indxOpt);

fopt = 2*1e-2*MW1*MT1 + 25*std(T, 1);

Amagnet = 2*MT1*MW1*1e-6;

%%
figure(1)
clf
hold on
JMAG_GA_torque = cell2mat(table2cell(readtable("data/JMAG_GA_torque.txt")));
angles = JMAG_GA_torque(:,1);
torque = JMAG_GA_torque(:,2);

plot(angles, torque);
grid on
xlabel("Rotation angle")
ylabel("Torque (Nm/m)")
set(gca, "FontSize", 12)

scaledStdT = 25*std(torque, 1);
StdT = std(torque,1)/4;

fopt1 = 2*1e-2*MW1*MT1 + StdT*1e2;

% str = "";
% for i = 1:29
%     str = str + "(T" + string(i) + "-" + "Tave)*(T" + string(i) + "-" + "Tave)+";
% end
