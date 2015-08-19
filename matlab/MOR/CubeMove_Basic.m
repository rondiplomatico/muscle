clear classes;

%% Init
% geo = 1; % Small
geo = 2; % Large
numsamples = 8;

% Processing
c = CubeMove('GeoNr',geo);
m = models.muscle.Model(c);
% Chosen so that num_dd*1.3 <= numsamples*T/dt
% to have 1.5 times more trajectory data than reducable dimensions.
m.dt = numsamples*m.T/(m.System.NumDerivativeDofs*1.5);

%m.EnableTrajectoryCaching = true;
%m.ComputeParallel = true;
%m.Data.useFileTrajectoryData;

m.Sampler = sampling.ManualSampler(linspace(.1,2,numsamples));
% This causes the 1st param to be used as training param
m.TrainingParams = 1;

d = approx.DEIM(m.System);
d.MaxOrder = m.System.NumDerivativeDofs;
%m.Approx = d;

%% Crunch
m.offlineGenerations;

%% build reduced
r = m.buildReducedModel(350);

%%
mu = m.DefaultMu;
mu(3) = .15; %new param!
[t,y] = m.simulate(mu,1);
[t,yr] = r.simulate(mu,1);
norm(y-yr)/norm(y)

%% DEIM stuff
%% Crunch
m.Approx = approx.DEIM(m.System);
m.off4_genApproximationTrainData;
m.Approx.MaxOrder = m.System.NumDerivativeDofs;
m.off5_computeApproximation;

%% check
r = m.buildReducedModel(350);
r.System.f.Order = 
[t,yr] = r.simulate(mu,1);

%% Plot
pm = PlotManager(false,2,2);
ma = ModelAnalyzer(r);
ma.plotReductionOverview(pm);





