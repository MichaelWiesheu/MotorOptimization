close all 
clear all 
clc
restoredefaultpath

addpath(genpath("../../"))
%%

angles = 0:1:29;
Tinit = load("../opt_param/plots/final/OptimizationHistory.mat").OptimizationHistory(1).Responses(end).Value *4;
Tparam = load("../opt_param/plots/final/OptimizationHistory.mat").OptimizationHistory(end).Responses(end).Value *4;

Tshape = load("../opt_shape/plots/final/OptimizationHistory.mat").OptimizationHistory(end).Responses(end).Value *4;
Tparamthenshape = load("../opt_param_then_shape/plots/final/OptimizationHistory.mat").OptimizationHistory(end).Responses(end).Value *4;
Tparamandshape = load("../opt_param_shape/plots/final/OptimizationHistory.mat").OptimizationHistory(end).Responses(end).Value *4;


Ainit = load("../opt_param/plots/final/OptimizationHistory.mat").OptimizationHistory(1).Responses(1).Value;
Aparam = load("../opt_param/plots/final/OptimizationHistory.mat").OptimizationHistory(end).Responses(1).Value;
Ashape = load("../opt_shape/plots/final/OptimizationHistory.mat").OptimizationHistory(end).Responses(1).Value;
Aparamthenshape = load("../opt_param_then_shape/plots/final/OptimizationHistory.mat").OptimizationHistory(end).Responses(1).Value;
Aparamandshape = load("../opt_param_shape/plots/final/OptimizationHistory.mat").OptimizationHistory(end).Responses(1).Value;


%%

allangles = load("data/Torque_init.mat").angles;
allTinit = load("data/Torque_init.mat").T*4;
allTparam = load("data/Torque_param.mat").T*4;
allTshape = load("data/Torque_shape.mat").T*4;
allTparamthenshape = load("data/Torque_param_then_shape.mat").T*4;
allTparamandshape = load("data/Torque_param_shape.mat").T*4;

lw = 1.5;
ms = 50;
figure(1)
clf
hold on
grid on

plot(allangles, allTinit, LineWidth=lw, Color="black");
plot(allangles, allTparam, LineWidth=lw, Color=TUDa_getColor('1a'));
plot(allangles, allTshape, LineWidth=lw, Color=TUDa_getColor('3a'));
plot(allangles, allTparamthenshape, LineWidth=lw, Color=TUDa_getColor('7a'));
plot(allangles, allTparamandshape, LineWidth=lw, Color=TUDa_getColor('9a'));

scatter(angles, Tparam, "filled", SizeData=ms, MarkerFaceColor=TUDa_getColor('1a'), Marker="o");
scatter(angles, Tshape, "filled", SizeData=ms, MarkerFaceColor=TUDa_getColor('3a'), Marker="s");
scatter(angles, Tparamthenshape, "filled", SizeData=ms, MarkerFaceColor=TUDa_getColor('7a'), Marker="^");
scatter(angles, Tparamandshape, "filled", SizeData=ms, MarkerFaceColor=TUDa_getColor('9a'), Marker="diamond");


ylim([1,2.7])
xlim([0,30])
ylabel("Torque (Nm)")
xlabel("Rotation angle (^\circ)")

legend(["Initial", "Parameter optimized", "Shape optimized", "Parameter then shape optimized", "Parameter and shape optimized"], Location="southeast")
set(gca,"FontSize", 18)


%% export to csv
allData(:,1) = allangles;
allData(:,2) = allTinit;
allData(:,3) = allTparam;
allData(:,4) = allTshape;
allData(:,5) = allTparamthenshape;
allData(:,6) = allTparamandshape;
writematrix(allData, "allData.csv")

optData(:,1) = angles;
optData(:,2) = Tparam;
optData(:,3) = Tshape;
optData(:,4) = Tparamthenshape;
optData(:,5) = Tparamandshape;
writematrix(optData, "optData.csv")



%%
disp(1-Aparam/Ainit)
disp(1-Ashape/Ainit)
disp(1-Aparamthenshape/Ainit)
disp(1-Aparamandshape/Ainit)
%%
disp(1-std(Tparam,1)/std(Tinit,1))
disp(1-std(Tshape,1)/std(Tinit,1))
disp(1-std(Tparamthenshape,1)/std(Tinit,1))
disp(1-std(Tparamandshape,1)/std(Tinit,1))