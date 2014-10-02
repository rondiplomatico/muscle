clear classes;
m = fullmuscle.Model(fullmuscle.CPull);
f = m.System.f;

mu0 = [1;0;0;0];
mu = [1;0;1;0]; 

%%
x0 = m.System.x0.evaluate(mu);
x00 = m.System.x0.evaluate(mu0);

m.System.prepareSimulation(mu,1);
%%
J = f.getStateJacobian(x0,0);
[JFD, dx] = f.getStateJacobianFD(x0,0);

%%
xin = 56;
fxout = 181;

dy = zeros(187,1);
h = 1e-2; %dx(xin);
dy(xin) = h;

fx = f.evaluate(x0,0);
fxh = f.evaluate(x0+dy,0);
dfx = (fxh-fx) / h;

%%
J(fxout,xin)
JFD(fxout,xin)

%%
fx(fxout)
fxh(fxout)
fxh(fxout)-fx(fxout)
dfx(fxout)
