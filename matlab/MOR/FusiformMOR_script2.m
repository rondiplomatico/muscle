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
%
%%
% build model
modconf = FusiformMORexample;
geo = modconf.PosFE.Geometry;
model = muscle.Model(modconf);
%
%%
% vary tolerances
model.ODESolver.RelTol = 0.5;
model.ODESolver.AbsTol = 0.5;
%
%%
% create parameter grid
s = sampling.GridSampler;
model.Sampler = s;
%
model.Data.TrajectoryData.UniformTrajectories = false;
model.Data.TrajectoryFxiData.UniformTrajectories = false;
%
model.off1_createParamSamples;
model.TrainingInputs = 1;
model.off2_genTrainingData;
model.save;
%
%%
% do offline phase 3 
%
% first simple
%model.SpaceReducer.MaxSubspaceSize = 1000;
%model.off3_computeReducedSpace;
%model.save;
%
% then test including fxi data
%IncludeTrajectoryFxiData = true;

