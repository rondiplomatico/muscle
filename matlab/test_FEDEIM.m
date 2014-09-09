clear classes;
mc = Cube12;
% mc = muscle.DebugConfig(8);
g = mc.PosFE.Geometry;
pg = mc.PressFE.Geometry;
m = muscle.Model(mc);
s = m.System;
f = s.f;

mu = m.DefaultMu;
[t,y,ct,x] = m.simulate(mu);

off_v_glob = 3*g.NumNodes;

nt = length(t);
fx = zeros(s.num_uvp_dof,nt);
fx_ass = fx;
fxu = zeros(f.num_uvp_dof_unass,nt);
for k = 1:nt
    f.ComputeUnassembled = true;
    fxu(:,k) = f.evaluate(x(:,k),t(k));
    f.ComputeUnassembled = false;
    fx(:,k) = f.evaluate(x(:,k),t(k));
end
% du part (no assembly)
fx_ass(1:s.num_u_dof,:) = fxu(1:s.num_u_dof,:);
% dvw part (with assembly)
fx_ass(s.num_u_dof + (1:s.num_v_dof+s.num_p_dof),:) = s.f.Sigma * fxu(s.num_u_dof+1:end,:);
diff = fx_ass - fx;
Norm.L2(diff)

% ass = sys.f.dvw_unass_elem_assoc;
% 
% 
% smc = muscle.SubMeshModelConfig(mc, [1 2]);
% sm = muscle.Model(smc);
% sm.plotGeometrySetup;