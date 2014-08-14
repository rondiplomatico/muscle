% script for the FusiformMORexample - start offline phase of MOR
%
% use GridSampler to create different parameter combinations
% do off1 and off2 (in BaseFullModel.m)
%
%%
clear all;
close all;
clear classes;
clc;
% build model
modconf = FusiformMORexample;
geo = modconf.PosFE.Geometry;
model = muscle.Model(modconf);
%
% vary tolerances
%%
model.ODESolver.RelTol = 1e-2;
model.ODESolver.AbsTol = 1e-2;
%%
% create parameter grid
% 
s = sampling.GridSampler;
model.Sampler = s;

model.off1_createParamSamples;
%model.Data.ParamSamples = model.Data.ParamSamples(:,1:3);
model.TrainingInputs = 1;
model.off2_genTrainingData;
model.save;

%save(fullfile(FusiformMORexample.OutputDir,'trajectoryData'),'model');
