clear classes;

outdir = '/home/dwirtz/software/muscle/data/';
pm = PlotManager;
pm.ExportDPI = 200;

%% %% Sparsity pattern (beispiel)
mc = fullmuscle.CPull(3,0:4:20,0:5:10);
m = fullmuscle.Model(mc);
ax = pm.nextPlot('sparsity_pattern','Overall system connectivity (Jacobian sparsity pattern)');
spy(m.System.f.JSparsityPattern);
set(get(gca,'Children'),'Parent',ax);
pm.savePlots(outdir,'Format','jpg','Close',true);

%% Motoneuron-Frequenz
load /home/dwirtz/software/muscle/data/oliver/motoneuron_freq.mat;
ax = pm.nextPlot('moto_freq','Motoneuron frequencies for different external signals and fibre types','Fibre type','External signal');
surf(ax,FT,MC,HZ,'FaceColor','interp','EdgeColor','interp');
view([-32 36]);
pm.savePlots(outdir,'Format','jpg','Close',true);

%%
base = '/home/dwirtz/software/muscle/data/oliver';
load('/home/dwirtz/software/muscle/data/cpull_v3_1000ms_50traj.mat');
ps = m.Data.ParamSamples;
[t,y] = m.simulate(ps(:,12),2);
m.plot(t,y,'Moto',false,'Freq',false,'Vid',fullfile(base,'cpull_extsig3.avi'));

[t,y] = m.simulate(ps(:,7),2);
m.plot(t,y,'Moto',false,'Freq',false,'Vid',fullfile(base,'cpull_balancedforce.avi'));

[t,y] = m.simulate(ps(:,2),1);
m.plot(t,y);

%% Fibre demo
% mc = fullmuscle.CFibreDemo;
% m = fullmuscle.Model(mc);
% [t,y] = m.simulate([1;0;0;2.5],2);
% [t,y] = m.simulate([1;0;0;4],2);
% [t,y] = m.simulate([1;0;0;9],1);
% [t,y] = m.simulate([1;0;0;9],2);
% [t,y] = m.simulate([1;0;0;9],3);
% m.plot(t,y);
% m.save;
% save /home/dwirtz/software/muscle/data/cfibredemo m;

