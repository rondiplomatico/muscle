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
load('/home/mordhorst/data/musclemodel__id-2482168/musclemodel_.mat');
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

save(fullfile(LongMORexample.OutputDir,'script_plotTraj_out'),'o');

% ggf unvollst. löschen
% incomplete = find(m.TrajectoriesCompleted(3,:) == 0);
% o(incomplete,:) = [];

%
%%
pm = PlotManager;
pm.LeaveOpen = true;
pm.FilePrefix = 'ypos_of_precomp_traj';
h = pm.nextPlot('plottag','Visualisation of  precomputed trajectories',...
            'time [ms]','y-position of right end');
% ALL
plot(h,t,o);
%
pm.done;
%
%%
% save
pm.SaveFormats = {'jpg','eps'};
pm.ExportDPI = 300;
pm.UseFileTypeFolders = false;
pm.NoTitlesOnSave = true;
pm.savePlots(LongMORexample.OutputDir);
%
%%