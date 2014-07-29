% script for testing the FusiformMORexample
clear all;
close all;
clear classes;
clc;
% build model
modconf = FusiformMORexample;
geo = modconf.PosFE.Geometry;
model = muscle.Model(modconf);
% simulate/calculate solution
[t,y] = model.simulate;
% plot solution
model.plot(t(1:10:end),y(:,1:10:end));
% get some output of interest
o = modconf.getOutputOfInterest(model,t,y);