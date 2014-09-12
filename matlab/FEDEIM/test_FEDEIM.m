% clear classes;
% muscle.Dynamics.test_UnassembledEvaluation;

%% Init
mc = Config_FEDEIM;
g = mc.PosFE.Geometry;
pg = mc.PressFE.Geometry;
m = muscle.Model(mc);
s = m.System;
f = s.f;

m.Sampler = sampling.ManualSampler;
m.Sampler.Samples = [1 1 5 5; 30 60 30 60];
m.TrainingParams = [1 2];
m.off1_createParamSamples;

%% Training data
m.off2_genTrainingData;
m.save;

%% FX-Evals
%ftd_fx_unass = data.FileTrajectoryData(fullfile(m.Data.DataDirectory,'fxtraj_unass'));
ftd_fx_unass = data.MemoryTrajectoryData;
f.ComputeUnassembled = true;
t = m.scaledTimes;
ntraj = m.Data.TrajectoryData.getNumTrajectories;
pi = ProcessIndicator('Computing unassembled evaluations for %d trajectories',ntraj,false,ntraj);
for k = 1:ntraj
    [x,mu] = m.Data.TrajectoryData.getTrajectoryNr(k);
    f.prepareSimulation(mu);
    tic;
    fx_u = f.evaluateMulti(x,t,mu);
    ct = toc;
    ftd_fx_unass.addTrajectory(fx_u,mu,[],ct);
    pi.step;
end
pi.stop;

%% Init 2
fxu = reshape(ftd_fx_unass.TrajectoryData,size(fxu,1),[]);
fx = zeros(s.num_uvp_dof,size(fxu,2));

% assemble full fx
fx(1:s.num_u_dof,:) = fxu(1:s.num_u_dof,:);
% dvw part (with assembly)
fx(s.num_u_dof + (1:s.num_v_dof+s.num_p_dof),:) = s.f.Sigma * fxu(s.num_u_dof+1:end,:);

save test;

d = general.DEIM;
fd = FEDEIM;
% 
% % %% Positions u
% % dof_u = fx(1:s.num_u_dof,:);
% % [U_dofu,S_dofu,V_dofu] = svd(dof_u,'econ');
% % pts_u = d.getInterpolationPoints(U_dofu);
% % err_u = fd.getInterpolErrors(U_dofu, pts_u, dof_u);
 
%% Classic DEIM
dof_vp = fx(s.num_u_dof+1:end,:);
[U,S,V] = svd(dof_vp,'econ');
num_DEIM = sum(diag(S) > S(1)*eps);
pts_DEIM = d.getInterpolationPoints(U(:,1:num_DEIM));
err_DEIM = fd.getInterpolErrors(U, pts_DEIM, dof_vp);
req_elems_DEIM = zeros(1,length(pts_DEIM));
for k = 1:length(pts_DEIM)
    req_elems_DEIM(k) = sum(sum(s.f.Sigma(pts_DEIM(1:k),:)*idx_elems') ~= 0);
end

%% UDEIM
dof_vp_unass = fxu(s.num_u_dof+1:end,:);
[Uu,Su,Vu] = svd(dof_vp_unass,'econ');
num_UDEIM = sum(diag(Su) > Su(1)*eps);
pts_UDEIM = d.getInterpolationPoints(Uu(:,1:num_UDEIM));
err_UDEIM = fd.getInterpolErrors(Uu, pts_UDEIM, dof_vp_unass);
% err_UDEIM_to_DEIM = fd.getInterpolErrors(Uu, pts_UDEIM, dof_vp_unass, s.f.Sigma);
req_elems_UDEIM = zeros(1,length(pts_UDEIM));
for k = 1:length(pts_UDEIM)
    req_elems_UDEIM(k) = sum(sum(s.f.Sigma(:,pts_UDEIM(1:k))*idx_elems(:,pts_UDEIM(1:k))') ~= 0);
end

%% FE-DEIM
% u = Uu(:,1:num_UDEIM);
% % u = Uu;
% % u = bsxfun(@times,Uu,1./max(Uu,[],2));
% % u = dof_vp_unass;
% % o = general.Orthonormalizer;
% % u = o.orthonormalize(u);
% [pts_FEDEIM, elems, basis_idx] = fd.getInterpolationElements(u, f.idx_vp_dof_unass_elems);
% err_FEDEIM = fd.getInterpolErrors(u, pts_FEDEIM, dof_vp_unass);
% % err_FEDEIM_to_DEIM = fd.getInterpolErrors(Uu, pts_FEDEIM, dof_vp_unass, s.f.Sigma);
% req_elems_FEDEIM = zeros(1,length(pts_FEDEIM));
% for k = 1:length(pts_FEDEIM)
%     req_elems_FEDEIM(k) = sum(sum(s.f.Sigma(:,pts_FEDEIM(1:k))*idx_elems(:,pts_FEDEIM(1:k))') ~= 0);
% end

%% FE-DEIM kNN
% u = Uu;
u = Uu(:,1:num_UDEIM);
allk = [1:10 12:2:20 25:30:size(u,1)];
if allk(end) ~= size(u,1)
    allk(end+1) = size(u,1);
end
nk = length(allk);
pts_kNN = cell(1,nk);
v_kNN = pts_kNN;
elems_kNN = pts_kNN;
errs_kNN = pts_kNN;
req_elems_kNN = pts_kNN;
% pi = ProcessIndicator('Crunching numbers',length(allk),false,length(allk));
for k = 1:nk
    fprintf('Starting k-NN with k=%d\n',k);
    [pts_loc, v_kNN{k}, elems_kNN{k}] = ...
        fd.getInterpolationPointskNN(u, allk(k), f.idx_vp_dof_unass_elems);
    pts_kNN{k} = pts_loc;
    errs_kNN{k} = fd.getInterpolErrors(u, pts_loc, dof_vp_unass, s.f.Sigma);
    req_elems = zeros(1,length(pts_loc));
    for idx = 1:length(pts_loc)
        req_elems(idx) = sum(sum(s.f.Sigma(:,pts_loc(1:idx))*idx_elems(:,pts_loc(1:idx))') ~= 0);
    end
    req_elems_kNN{k} = req_elems;
%     pi.step;
end
pi.stop;

%% FE-DEIM Tests
% u = Uu(:,1:num_UDEIM);
% u = Uu;
% u = bsxfun(@times,Uu,1./max(Uu,[],2));
% u = dof_vp_unass;
% o = general.Orthonormalizer;
% u_orth = o.orthonormalize(u);

% pts_FEDEIM = pts;
% err_FEDEIM = fd.getInterpolErrors(Uu, pts_FEDEIM, dof_vp_unass);
% % err_FEDEIM_to_DEIM = fd.getInterpolErrors(Uu, pts_FEDEIM, dof_vp_unass, s.f.Sigma);
% req_elems_FEDEIM = zeros(1,length(pts_FEDEIM));
% for k = 1:length(pts_FEDEIM)
%     req_elems_FEDEIM(k) = sum(sum(s.f.Sigma(:,pts_FEDEIM(1:k))*idx_elems(:,pts_FEDEIM(1:k))') ~= 0);
% end

%% Plots
pm = PlotManager(false,3,2);
pm.AutoTickMarks = false;
% fd.plot(pm,'u',S_dofu,err_u);
fd.plot(pm,'DEIM',S,err_DEIM,req_elems_DEIM);
fd.plot(pm,'UDEIM',Su,err_UDEIM,req_elems_UDEIM);
% fd.plot(pm,'UDEIM_DEIM',Su,err_UDEIM_to_DEIM);
for k = 1:length(allk)
    fd.plot(pm,sprintf('UDEIM-%dNN',allk(k)),Su,errs_kNN{k},req_elems_kNN{k});
end
% fd.plot(pm,'FEDEIM',Su,err_FEDEIM,req_elems_FEDEIM);
% fd.plot(pm,'FEDEIM',Su,err_FEDEIM_to_DEIM);
pm.done;

%% kNN - Evaluation
allk = allk(1:end-1);
logerrs = 1:.5:13; %10e-n
nk = length(allk);
ne = length(logerrs);
num_elems = zeros(nk,ne);
for k = 1:nk
    err = errs_kNN{k};
    for eidx = 1:ne
        pos = find(err(3,:) < 10^(-logerrs(eidx)),1,'first');
        if isempty(pos)
            num_elems(k,eidx) = 0;
        else
            num_elems(k,eidx) = req_elems_kNN{k}(pos);
        end
    end
end
ax = gca;
[X,Y] = meshgrid(10.^-logerrs,allk);
surf(X,Y,num_elems,'FaceColor','interp','EdgeColor','interp','Parent',ax);
xlabel('Max error below');
set(ax,'XScale','log');
ylabel('k');










