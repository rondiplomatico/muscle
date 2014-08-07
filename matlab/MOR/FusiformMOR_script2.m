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
file = fullfile(FusiformMORexample.OutputDir,'FusiformMOR_secondtest.mat');
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
% create parameter grid
% 
s = sampling.GridSampler;
model.Sampler = s;

model.off1_createParamSamples;

model.TrainingInputs = 1;
model.off2_genTrainingData;

