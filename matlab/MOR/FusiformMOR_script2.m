%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% script for the FusiformMORexample - start offline phase of MOR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% clear all;
% close all;
% clear classes;
% clc;
%
%%
% build model
%
% modconf = FusiformMORexample;
% geo = modconf.PosFE.Geometry;
% model = muscle.Model(modconf);
%
%%
% vary tolerances
%
% model.ODESolver.RelTol = 1e-3;
% model.ODESolver.AbsTol = 1e-3;
%
%%
% do offline phase 1&2
%
% create parameter grid
%
% s = sampling.GridSampler;
% model.Sampler = s;
% %
% model.Data.TrajectoryData.UniformTrajectories = false;
% model.Data.TrajectoryFxiData.UniformTrajectories = false;
% %
% model.off1_createParamSamples;
% %
% % compute trajectories and trajectories_fxi
% %
% model.TrainingInputs = 1;
% model.ComputeParallel = true;
% model.off2_genTrainingData;
% model.save;
%
%%
% do offline phase 3 
%
% model.SpaceReducer.MaxSubspaceSize = 1000;
% model.SpaceReducer.MinRelImprovement = 1e-7;
% model.off3_computeReducedSpace;
% model.save;
%
% test including fxi data
%IncludeTrajectoryFxiData = true;
%
%%
% do offline phase 4
%
sel = data.selection.LinspaceSelector;
sel.EnsureUniqueData = true;
sel.Size = 40000;
a = approx.DEIM(model.System);
a.TrainDataSelector = sel;
model.Approx = a;
model.off4_genApproximationTrainData;
model.save;
%
%%
% offline phase 5
%
a.MaxOrder = 1000;
model.off5_computeApproximation;
model.save;
%
%%
%
