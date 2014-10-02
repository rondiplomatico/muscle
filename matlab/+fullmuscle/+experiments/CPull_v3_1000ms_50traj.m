clear classes;
m = fullmuscle.Model(fullmuscle.CPull(3));

%% Example data setup
m.Sampler = sampling.ManualSampler;
m.TrainingParams = [3 4];
m.Sampler.Samples = Utils.createCombinations([1 5 10 15 20],[0 1 3 5 9]);
m.TrainingInputs = [1 2];
m.Data.useFileTrajectoryData;
m.Data.TrajectoryData.UniformTrajectories = false;

%% Offline
m.off1_createParamSamples;
matlabpool open;
m.ComputeParallel = true;
m.off2_genTrainingData;
matlabpool close;

m.save;
save /home/dwirtz/software/muscle/data/cpull_v3_1000ms_50traj m;