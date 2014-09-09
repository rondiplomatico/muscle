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
for k = 1:nt
    f.ComputeUnassembled = true;
    fxu(:,k) = f.evaluate(x(:,k),t(k));
    f.ComputeUnassembled = false;
    fx(sys.dof_idx_global,k) = f.evaluate(x(:,k),t(k));
end
% du part (no assembly)
ndirdispl = length(sys.dof_idx_displ);
fx_ass(sys.dof_idx_displ,:) = fxu(1:ndirdispl,:);
% dvw part (with assembly)
idx = sys.dof_idx_global;
idx(1:ndirdispl) = [];
fx_ass(idx,:) = sys.f.Sigma * fxu(off_v_glob+1:end,:);
diff = fx_ass - fx;
Norm.L2(diff)

% ass = sys.f.dvw_unass_elem_assoc;
% 
% 
% smc = muscle.SubMeshModelConfig(mc, [1 2]);
% sm = muscle.Model(smc);
% sm.plotGeometrySetup;