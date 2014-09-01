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
% load('/home/mordhorst/data/musclemodel__id-2129151/musclemodel_.mat');%0.5 LEAD
% load('/data/homes/mordhorst/data/musclemodel__id-2129151/musclemodel_.mat');%0.5 local
% load('/home/mordhorst/data/musclemodel__id-2140558/musclemodel_.mat');%1e-3
load('/home/mordhorst/data/musclemodel__id-2112992/musclemodel_.mat');%0.1
%
pi = ProcessIndicator('get trajectory data',1,false,1);
j = 1;
for k = 1:125
%     if ~any(k == [1 11 16 21 31 36 41 46 81 86])
%     if ~any(k == [1 6 11 16 21 46 51 56 61 66 71 76 78 81 86 91 96 101])
    if ~any(k == [6 11 16 21 28 31 36 41 46 56 61 71 91])
        [t,y] = m.simulate(m.Data.ParamSamples(:,k));
        o(j,:) = m.Config.getOutputOfInterest(m,t,y);
%         o(k,:) = m.Config.getOutputOfInterest(m,t,y);
        j = j+1;
        pi.step(1/115);
    end
end
pi.stop;
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
pm.SaveFormats = {'jpg'};
pm.ExportDPI = 300;
pm.UseFileTypeFolders = false;
pm.NoTitlesOnSave = true;
%pm.savePlots('/home/mordhorst/data/MORresults');%LEAD
pm.savePlots('/data/homes/mordhorst/data/MORresults');%local
%
%%
