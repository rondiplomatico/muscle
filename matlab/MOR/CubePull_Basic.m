clear classes;

%% Init
% geo = 1; % Small
geo = 2; % Large
numsamples = 5;

% Processing
c = models.muscle.examples.CubePull('GeoNr',geo);
m = models.muscle.Model(c);
% Chosen so that num_dd*1.3 <= numsamples*T/dt
% to have 1.5 times more trajectory data than reducable dimensions.
m.dt = numsamples*m.T/(m.System.NumDerivativeDofs*1.5);

%m.EnableTrajectoryCaching = true;
%m.ComputeParallel = true;
%m.Data.useFileTrajectoryData;

forces = linspace(-.1,.5,numsamples);
m.Sampler = sampling.ManualSampler(forces);
% This causes the 3rd param to be used as training param
m.TrainingParams = 3;
m.TrainingInputs = 1;

d = approx.DEIM(m.System);
d.MaxOrder = m.System.NumDerivativeDofs;
%m.Approx = d;

%% Crunch
m.offlineGenerations;

%% Test algebraic constraint violations for subspace sizes
mu = m.DefaultMu;
% x0 = m.System.getX0(mu);
ad = m.System.NumAlgebraicDofs;
% m.System.prepareSimulation(mu,rm.DefaultInput);
res = zeros(1, m.System.NumStateDofs);
for ridx = 1:10:m.System.NumStateDofs
%     rx0 = rm.System.R*rm.System.R'*x0;
%     yp0 = rode(0,x0);
%     constr_violation = yp0(end-ad+1:end);
%     res(r) = Norm.L2(constr_violation);
    r = m.buildReducedModel(ridx);
    x0 = r.System.getX0(mu);
    r.System.prepareSimulation(mu,r.DefaultInput);
    yp0 = r.System.ODEFun(0,x0);
    constr_violation = yp0(end-ad+1:end);
    res(ridx) = Norm.L2(constr_violation);
end
semilogy(res);

%% build reduced
% r = m.buildReducedModel(368);
r = m.buildReducedModel(200);

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
r = m.buildReducedModel(200);
[t,yr] = r.simulate(mu,1);

%% Plot
pm = PlotManager(false,2,2);
ma = ModelAnalyzer(r);
ma.plotReductionOverview(pm);





