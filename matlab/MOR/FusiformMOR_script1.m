% script for testing the FusiformMORexample
clear all;
close all;
clear classes;
clc;
% build model
file = fullfile(FusiformMORexample.OutputDir,'FusiformMOR_firsttests.mat');
saveit = true;
if exist(file,'file') == 2
    load(file);
    saveit = false;
else
    modconf = FusiformMORexample;
    geo = modconf.PosFE.Geometry;
    model = muscle.Model(modconf);
end
%
%%
pi1 = ProcessIndicator('viscosity test',1,false,1);
%disp('start viscosity tests');
o_visc = zeros(100,size(model.Times,2));
j=1;
for k = 0.01:0.1:10
    %k
    % simulate/calculate solution
    [t,y] = model.simulate([k;50;-10]);
    % plot solution
    % model.plot(t(1:10:end),y(:,1:10:end));
    % get some output of interest
    o_visc(j,:) = modconf.getOutputOfInterest(model,t,y);
    j=j+1;
    pi1.step(1/100);
end
pi1.stop;
%disp('simulation done')
%
%%
% plot data
figure()
plot(t,o_visc);
xlabel('time [ms]');
ylabel('mean position of nodes at right end in y-direction');
% legend('viscosity=1','viscosity=2','viscosity=3','viscosity=4','viscosity=5',...
%         'viscosity=6','viscosity=7','viscosity=8','viscosity=9','viscosity=10')
title('apply force for 10ms up to -10, after 20ms activate for 50ms from 0>1, different viscosities in [0.01,10]')    
%
%%
pi2 = ProcessIndicator('activation test',1,false,1);
%disp('start activation tests');
o_alpha = zeros(20,size(model.Times,2));
j=1;
for k = 10:10:200
    %k
    % simulate/calculate solution
    [t,y] = model.simulate([1;k;-5]);
    % plot solution
    % model.plot(t(1:10:end),y(:,1:10:end));
    % get some output of interest
    o_alpha(j,:) = modconf.getOutputOfInterest(model,t,y);
    j=j+1;
    pi2.step(1/20);
end
pi2.stop;
%disp('simulation done')
%
%%
% plot data
figure()
plot(t,o_alpha(1:2:end,:));
xlabel('time [ms]');
ylabel('mean position of nodes at right end in y-direction');
legend('alpha-ramp=10ms','alpha-ramp=30ms','alpha-ramp=50ms','alpha-ramp=70ms',...
        'alpha-ramp=90ms','alpha-ramp=110ms','alpha-ramp=130ms','alpha-ramp=150ms')
title('apply force for 10ms up to -5, after 20ms activate for value legend ms from 0>1, viscosity=1')    
%
%%
pi3 = ProcessIndicator('force test',1,false,1);
%disp('start force tests');
o_force = zeros(16,size(model.Times,2));
j=1;
for k = 0:10:150
    %k
    % simulate/calculate solution
    [t,y] = model.simulate([1;50;-k]);
    % plot solution
    % model.plot(t(1:10:end),y(:,1:10:end));
    % get some output of interest
    o_force(j,:) = modconf.getOutputOfInterest(model,t,y);
    j=j+1;
    pi3.step(1/16);
end
pi3.stop;
%disp('simulation done')
%
%%
% plot data
figure()
plot(t,o_force);
xlabel('time [ms]');
ylabel('mean position of nodes at right end in y-direction');
legend('force=0','force=-10','force=-20','force=-30','force=-40',...
        'force=-50','force=-60','force=-70','force=-80','force=-90','force=-100',...
        'force=-110','force=-120','force=-130','force=-140','force=-150')
title('apply force for 10ms up to value legend, after 20ms activate for 50ms from 0>1, viscosity=1')
%
%%
pi4 = ProcessIndicator('extremes test',1,false,1);
%disp('start extremes tests');
%o = zeros(8,size(model.Times,2));
% viscosity range [0.01,10]
mu1_min = 1e-3;
mu1_max = 5;%10;
% activation duration range [10ms, 200ms]
mu2_min = 1;
mu2_max = 200;
% force range [0,-150]
mu3_min = 0;
mu3_max = -1e3;
%
extremes_tests(:,1) = [mu1_min; mu2_min; mu3_min];
extremes_tests(:,2) = [mu1_min; mu2_min; mu3_max];
extremes_tests(:,3) = [mu1_min; mu2_max; mu3_min];
extremes_tests(:,4) = [mu1_min; mu2_max; mu3_max];
extremes_tests(:,5) = [mu1_max; mu2_min; mu3_min];
extremes_tests(:,6) = [mu1_max; mu2_min; mu3_max];
extremes_tests(:,7) = [mu1_max; mu2_max; mu3_min];
extremes_tests(:,8) = [mu1_max; mu2_max; mu3_max];
% test extremes
for k = 1:8
    %k
    % simulate/calculate solution
    [t,y] = model.simulate(extremes_tests(:,k));
    % get some output of interest
    o(k,:) = modconf.getOutputOfInterest(model,t,y);
    pi4.step(1/8);
end
pi4.stop;
%disp('simulation done')
%
%
%%
% plot data
figure()
plot(t,o);
xlabel('time [ms]');
ylabel('mean position of nodes at right end in y-direction');
legend('[min,min,min]','[min,min,max]','[min,max,min]','[min,max,max]',...
        '[max,min,min]','[max,min,max]','[max,max,min]','[max,max,max]')
title('extreme tests: viscosity mu_1, activation starting at 20ms with duration mu_2, force for 10ms up to mu_3')
%%
if saveit
    save(file,'model','modconf','t','y','o_visc','o_force','o_alpha','o');
end
    