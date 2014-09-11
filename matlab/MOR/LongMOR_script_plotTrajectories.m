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
% load('/home/mordhorst/data/musclemodel__id-2482168/musclemodel_.mat');%LEAD
load('/data/homes/mordhorst/data/musclemodel__id-2482168/musclemodel_.mat');%local
%
np = m.Data.SampleCount;
oy = -Inf(np,length(m.Times));
oz = -Inf(np,length(m.Times));
pi = ProcessIndicator('Reading %d trajectories',np,false,np);
for k = 1:np
    [t,y] = m.simulate(m.Data.ParamSamples(:,k));
    [oy(k,1:length(t)),oz(k,1:length(t))] = m.Config.getOutputOfInterest(m,t,y);
    pi.step;
end
pi.stop;
%
save(fullfile(LongMORexample.OutputDir,'script_plotTraj_out'),'oy','oz');
%
%%
pm = PlotManager;
pm.LeaveOpen = true;
pm.FilePrefix = 'pos_of_precomp_traj';
%
% z-position of middle node
h = pm.nextPlot('zpos','Visualisation of  precomputed trajectories',...
            'time [ms]','z-position of middle node');
plot(h,t,oz);
%
pm.done;
%
%%
pm = PlotManager;
pm.LeaveOpen = true;
pm.FilePrefix = 'pos_of_precomp_traj';
%
% y-position of node at 1/3rd y-direction
h = pm.nextPlot('ypos','Visualisation of  precomputed trajectories',...
            'time [ms]','y-position of node (0,-20,0)');
plot(h,t,oy);
%
pm.done;
%
%%
pm = PlotManager;
pm.LeaveOpen = true;
pm.FilePrefix = 'pos_of_precomp_traj';
%
% activations for viscosity 0.5 OR 5.1 OR 10
h = pm.nextPlot('zpos_visc0p5_act','Visualisation of  precomputed trajectories',...
            'time [ms]','z-position of middle node',{'activation in 20ms',...
            'activation in 24ms','activation in 28ms','activation in 34ms',...
            'activation in 41ms','activation in 49ms','activation in 58ms',...
            'activation in 69ms','activation in 83ms','activation in 100ms','Location','SouthEast'});
plot(h,t,oz(1:10:100,:));% visc 0.5 - 1, visc 5.1 - 8; visc 10 - 10
%
pm.done;
%
%%
pm = PlotManager;
pm.LeaveOpen = true;
pm.FilePrefix = 'pos_of_precomp_traj';
%
% normalised activations for viscosity 0.5 OR 5.1 OR 10
h = pm.nextPlot('zpos_visc0p5_act-normalised','Visualisation of  precomputed trajectories',...
            'activation','normalised z-position of middle node',{'activation in 20ms',...
            'activation in 24ms','activation in 28ms','activation in 34ms',...
            'activation in 41ms','activation in 49ms','activation in 58ms',...
            'activation in 69ms','activation in 83ms','activation in 100ms','Location','SouthEast'});
sel = 1:10:100;% visc 0.5 - 1, visc 5.1 - 8; visc 10 - 10
act_times = m.Data.ParamSamples(2,sel);
cc = hsv(10);
for i = 1:10
    axis([0 1 0 8]);
    plot(h,t/act_times(i),oz(sel(i),:),'color',cc(i,:));
    hold on;
end
hold off;
%
pm.done;
%
%%
pm = PlotManager;
pm.LeaveOpen = true;
pm.FilePrefix = 'pos_of_precomp_traj';
%
% plot over parameter space
%
paramspace_visc = reshape(m.Data.ParamSamples(1,:)',10,10);
paramspace_act = reshape(m.Data.ParamSamples(2,:)',10,10);
%
h = pm.nextPlot('ypos','Visualisation of  precomputed trajectories',...
    'time to full activation','viscosity');
zlabel(h,'z-position of middle node');
axis(h,[20 100 0.5 10 0 8]);
hold(h,'on');
%
for i = 1:101
    mesh_oz = reshape(oz(:,i),10,10);
    cla(h);
    %
    mesh(h,paramspace_act,paramspace_visc,mesh_oz);
    pause(.1);%drawnow;
end
%
pm.done;
%
%
%%
% plot or video over parameter space
video = false;
%
pm = PlotManager;
pm.LeaveOpen = true;
pm.FilePrefix = 'pos_of_precomp_traj';
%
if video
    avifile = fullfile(LongMORexample.OutputDir,'parameterspace.avi');
    vw = VideoWriter(avifile);
    vw.FrameRate = 25;
    vw.open;
end
%
h = pm.nextPlot('zpos',sprintf('z-position of middle node at t=%g',t(end)),...
    'time to full activation','viscosity');
zlabel(h,'z-position of middle node');
axis(h,[20 100 0.5 10 0 8]);
hold(h,'on');
%
paramspace_visc = reshape(m.Data.ParamSamples(1,:)',10,10);
paramspace_act = reshape(m.Data.ParamSamples(2,:)',10,10);
%
for ts = 1:length(t)
    mesh_oz = reshape(oz(:,ts),10,10);
    cla(h);
    %
    mesh(h,paramspace_act,paramspace_visc,mesh_oz);
    title(h,sprintf('z-position of middle node at t=%g',t(ts)));
    %
    if video
        vw.writeVideo(getframe(gcf));
    else
        pause(.2);
    end
end
%
if video
    vw.close;
end
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