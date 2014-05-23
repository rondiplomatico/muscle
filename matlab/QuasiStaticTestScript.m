m = muscle.Model(QuasiStaticTest);
m.T = 150;
m.dt = 1;
f = m.System.f;
f.alpha = 1;
% f.Pmax = 250;
f.rampFraction = .2;
m.System.Viscosity = 0;
os = m.ODESolver;
os.RelTol = .1;
os.AbsTol = .1;
[t,y] = m.simulate(1);
m.plot(t,y,'Velo',true);