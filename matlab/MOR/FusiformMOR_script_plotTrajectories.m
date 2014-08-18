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
%load('/home/mordhorst/data/musclemodel__id-2129151/musclemodel_.mat');%LEAD
load('/data/homes/mordhorst/data/musclemodel__id-2129151/musclemodel_.mat');%local
%
pi = ProcessIndicator('get trajectory data',1,false,1);
j = 1;
for k = 1:125
    if ~any(k == [1 11 16 21 31 36 41 46 81 86])
        [t,y] = m.simulate(m.Data.ParamSamples(:,k));
        o(j,:) = m.Config.getOutputOfInterest(m,t,y);
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
            'time [ms]','mean position of nodes at right end in y-direction');
plot(h,t,o);
pm.done;
%
%%
% save
pm.SaveFormats = {'jpg'};
pm.JPEGQuality = '95';
pm.UseFileTypeFolders = false;
%pm.savePlots('/home/mordhorst/data/MORresults');%LEAD
pm.savePlots('/data/homes/mordhorst/data/MORresults');%local
%
%%