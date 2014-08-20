%---------------------------------------------------------------------%
%
% script to visualise precomputed trajectories
%
%---------------------------------------------------------------------%
%%
clear all;
close all;
clear classes;
clc;
%
%%
% load('/home/mordhorst/data/musclemodel__id-2129151/musclemodel_.mat');%0.5
% load('/home/mordhorst/data/musclemodel__id-2140558/musclemodel_.mat');%1e-3
% load('/home/mordhorst/data/musclemodel__id-2112992/musclemodel_.mat');%0.1
load('/home/dwirtz/data/musclemodel__id-2101435/musclemodel_.mat');
%

np = m.Data.SampleCount;
o = -Inf(np,length(m.Times));
pi = ProcessIndicator('Reading %d trajectories',np,false,np);
for k = 1:np
    [t,y] = m.simulate(m.Data.ParamSamples(:,k));
    o(k,1:length(t)) = m.Config.getOutputOfInterest(m,t,y);
    pi.step;
end
pi.stop;

save(fullfile(FusiformMORexample.OutputDir,'script_plotTraj_out'),'o');

% ggf unvollst. löschen
% incomplete = find(m.TrajectoriesCompleted(3,:) == 0);
% o(incomplete,:) = [];

%
%%
pm = PlotManager;
pm.LeaveOpen = true;
pm.FilePrefix = 'ypos_of_precomp_traj';
h = pm.nextPlot('plottag','Visualisation of  precomputed trajectories (115)',...
            'time [ms]','y-position of right end');
% ALL
plot(h,t,o);
%
% vary ACTIVATION, fixed viscosity and all forces
%plot(h,t,o(2:5:25,:));% mu = [2.5;act;-1000]
%plot(h,t,o(27:5:50,:));% mu = [2.5;act;-750]
%plot(h,t,o(52:5:75,:));% mu = [2.5;act;-500]
%plot(h,t,o(77:5:100,:));% mu = [2.5;act;-250]
%plot(h,t,o(102:5:125,:));% mu = [2.5;act;0]
%
% vary FORCE, fixed viscosity and all activation times
%plot(h,t,o([3:25:125],:));% mu = [5;1;force]
%plot(h,t,o([8:25:125],:));% mu = [5;50.75;force]
%plot(h,t,o([13:25:125],:));% mu = [5;100.5;force]
%plot(h,t,o([18:25:125],:));% mu = [5;150.25;force]
%plot(h,t,o([23:25:125],:));% mu = [5;200;force]
%
pm.done;
%
%%
% save
pm.SaveFormats = {'jpg','eps'};
pm.ExportDPI = 300;
pm.UseFileTypeFolders = false;
pm.NoTitlesOnSave = true;
pm.savePlots(FusiformMORexample.OutputDir);
%
%%