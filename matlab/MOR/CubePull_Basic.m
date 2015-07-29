clear classes;

%% Init
m = muscle.Model(CubePull);
%m.EnableTrajectoryCaching = true;

%m.ComputeParallel = true;
%m.Data.useFileTrajectoryData;

forces = linspace(-.1,.5,2);
m.Sampler = sampling.ManualSampler(forces);
% This causes the 3rd param to be used as training param
m.TrainingParams = 3;
m.TrainingInputs = 1;

d = approx.DEIM(m.System);
d.MaxOrder = m.System.NumStateDofs;
%m.Approx = d;

%% Crunch
m.offlineGenerations;

%%
r = m.buildReducedModel(368);

%%
mu = m.DefaultMu;
mu(3) = .15; %new param!
[t,y] = m.simulate(mu,1);
[t,yr] = r.simulate(mu,1);
norm(y-yr)/norm(y)

%% Plot
pm = PlotManager(false,2,2);
ma = ModelAnalyzer(r);
ma.plotReductionOverview(pm);

%% Manual subspace comp
A = m.Data.TrajectoryData.toMemoryMatrix;
U = A(1:610,:);
V = A(611:1220,:);
[uu,us,uv] = svd(U,'econ');
[vu,vs,vv] = svd(V,'econ');

%% Plot manual
ax = pm.nextPlot('singvals_split',...
    'Singular values of same training data when split into u/v parts','singular value number','value');
semilogy(ax,1:610,diag(us),'r',1:610,diag(vs),'b')
pm.done;
pm.savePlots(pwd,'Format',{'jpg','pdf'});