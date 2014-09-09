clear classes;
mc = Cube12;
% mc = muscle.DebugConfig(8);
g = mc.PosFE.Geometry;
pg = mc.PressFE.Geometry;
m = muscle.Model(mc);
sys = m.System;
f = sys.f;

mu = m.DefaultMu;
[t,y,ct,x] = m.simulate(mu);

off_v_glob = 3*g.NumNodes;

nt = length(t);
fx = zeros(6*g.NumNodes+pg.NumNodes,nt);
fx_ass = fx;
dofs_unass = off_v_glob + 3*g.NumElements*g.DofsPerElement...
    + pg.NumElements*pg.DofsPerElement ...
    - length(f.bc_dir_idx_unass);
fxu = zeros(dofs_unass,nt);
bcidx = [off_v_glob+sys.bc_dir_displ_idx; sys.bc_dir_velo_idx];
for k = 1:nt
    f.ComputeUnassembled = true;
    fxu(:,k) = f.evaluate(x(:,k),t(k));
    f.ComputeUnassembled = false;
    fx(sys.dof_idx_global,k) = f.evaluate(x(:,k),t(k));
    %fx(bcidx,k) = f.LastBCResiduals;
end
fx_ass(sys.dof_idx_global,:) = sys.f.Sigma * fxu;
diff = fx_ass - fx;
Norm.L2(diff)

smc = muscle.SubMeshModelConfig(mc, [1 2]);
sm = muscle.Model(smc);
sm.plotGeometrySetup;