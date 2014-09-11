% % clear classes;
% % muscle.Dynamics.test_UnassembledEvaluation;
% 
% mc = Cube12;
% % mc = LongForceBC;
% % mc = EntireTA;
% % mc = muscle.DebugConfig(8);
% g = mc.PosFE.Geometry;
% pg = mc.PressFE.Geometry;
% m = muscle.Model(mc);
% s = m.System;
% f = s.f;
% 
% m.Sampler = sampling.ManualSampler;
% m.Sampler.Samples = [1 1 5 5; 10 30 10 30];
% m.TrainingParams = [1 2];
% m.off1_createParamSamples;
% 
% m.off2_genTrainingData;
% 
% %ftd_fx_unass = data.FileTrajectoryData(fullfile(m.Data.DataDirectory,'fxtraj_unass'));
% ftd_fx_unass = data.MemoryTrajectoryData;
% f.ComputeUnassembled = true;
% t = m.scaledTimes;
% ntraj = m.Data.TrajectoryData.getNumTrajectories;
% pi = ProcessIndicator('Computing unassembled evaluations for %d trajectories',ntraj,false,ntraj);
% for k = 1:ntraj
%     [x,mu] = m.Data.TrajectoryData.getTrajectoryNr(k);
%     f.prepareSimulation(mu);
%     tic;
%     fx_u = f.evaluateMulti(x,t,mu);
%     ct = toc;
%     ftd_fx_unass.addTrajectory(fx_u,mu,[],ct);
%     pi.step;
% end
% pi.stop;
% 
% fxu = reshape(ftd_fx_unass.TrajectoryData,size(fxu,1),[]);
% fx = zeros(s.num_uvp_dof,size(fxu,2));
% 
% % assemble full fx
% fx(1:s.num_u_dof,:) = fxu(1:s.num_u_dof,:);
% % dvw part (with assembly)
% fx(s.num_u_dof + (1:s.num_v_dof+s.num_p_dof),:) = s.f.Sigma * fxu(s.num_u_dof+1:end,:);
% 
% d = general.DEIM;
% fd = FEDEIM;
% 
% % %% Positions u
% % dof_u = fx(1:s.num_u_dof,:);
% % [U_dofu,S_dofu,V_dofu] = svd(dof_u,'econ');
% % pts_u = d.getInterpolationPoints(U_dofu);
% % err_u = fd.getInterpolErrors(U_dofu, pts_u, dof_u);
% 
% %% Classic DEIM
% dof_vp = fx(s.num_u_dof+1:end,:);
% [U,S,V] = svd(dof_vp,'econ');
% num_DEIM = sum(diag(S) > S(1)*eps);
% pts_DEIM = d.getInterpolationPoints(U(:,1:num_DEIM));
% err_DEIM = fd.getInterpolErrors(U, pts_DEIM, dof_vp);
% elems_DEIM = zeros(1,length(pts_DEIM));
% for k = 1:length(pts_DEIM)
%     elems_DEIM(k) = sum(sum(s.f.Sigma(pts_DEIM(1:k),:)*idx_elems') ~= 0);
% end
% 
% %% UDEIM
% dof_vp_unass = fxu(s.num_u_dof+1:end,:);
% [Uu,Su,Vu] = svd(dof_vp_unass,'econ');
% num_UDEIM = sum(diag(Su) > Su(1)*eps);
% pts_UDEIM = d.getInterpolationPoints(Uu(:,1:num_UDEIM));
% err_UDEIM = fd.getInterpolErrors(Uu, pts_UDEIM, dof_vp_unass);
% err_UDEIM_to_DEIM = fd.getInterpolErrors(Uu, pts_UDEIM, dof_vp_unass, s.f.Sigma);
% elems_UDEIM = zeros(1,length(pts_UDEIM));
% for k = 1:length(pts_UDEIM)
%     elems_UDEIM(k) = sum(sum(s.f.Sigma(:,pts_UDEIM(1:k))*idx_elems(:,pts_UDEIM(1:k))') ~= 0);
% end
% 
% %% Plots
% pm = PlotManager(false,3,2);
% pm.AutoTickMarks = false;
% % fd.plot(pm,'u',S_dofu,err_u);
% fd.plot(pm,'DEIM',S,err_DEIM);
% fd.plot(pm,'UDEIM',Su,err_UDEIM);
% fd.plot(pm,'UDEIM_DEIM',Su,err_UDEIM_to_DEIM);
% pm.done;

%% FE-DEIM
% elems_FEDEIM = d.getInterpolationElements(Uu,f.idx_vp_dof_unass_elems);
u = Uu(:,1:num_UDEIM);
% u = bsxfun(@times,Uu,1./max(Uu,[],2));
% u = dof_vp_unass;
% o = general.Orthonormalizer;
% u_orth = o.orthonormalize(u);
idx_elems = f.idx_vp_dof_unass_elems;

sigma = double(idx_elems);
num_elems = size(sigma,1);
n = size(u,1);
elems = zeros(1,num_elems);
maxerr = zeros(1,num_elems+1);
maxerr_res = maxerr;

% Init
resi_idx = 1:size(u,2);
resi_pos = true(1,size(u,2));
residual = u;
pts = [];
basis_idx = [];
basis = [];
% hlp = u;
v = [];
aborted = false;
for elem_idx = 1:num_elems
    % Determine element error measure
    no = sqrt(sigma*(residual.^2));
    [maxerr(elem_idx), elems(elem_idx)] = max(max(no,[],2));
    maxerr_res(elem_idx) = max(abs(residual(:)));
    % Get corresponding points
    pts_elem = find(idx_elems(elems(elem_idx),:));
    P_elem = sparse(pts_elem,1:length(pts_elem),ones(length(pts_elem),1),n,length(pts_elem));
%     
%     % Do a local DEIM point search on restricted set
%     [v(end+1), pt_loc] = max(abs(u(pts_elem,1)));
%     pts = [pts pts_elem(pt_loc)];%#ok
%     npts = length(pts);
%     P = sparse(pts,1:npts,ones(npts,1),n,npts);
%     v_loc_max = v(end);
%     for i = 2:length(pts_elem)
%         if v(end) < v_loc_max*sqrt(eps);
%             break;
%         end
%         c = (P'*u(:,1:(i-1))) \ (P'*u(:,i));
%         [v(end+1), pt_loc] = max(abs(u(pts_elem,i) - u(pts_elem,1:(i-1))*c));
%         pts = [pts pts_elem(pt_loc)];%#ok
%         npts = length(pts);
%         P = sparse(pts,1:npts,ones(npts,1),n,npts);
%     end
%     resi_pos(pts) = false;
%     c = (P'*u(:,1:npts)) \ (P'*u(:,resi_pos));
%     residual = u(:,resi_pos) - u(:,1:npts)*c;
    
    % Do a local EIM point search on restricted set
    [pts_local, basis_local, v] = fd.getInterpolationPointsMaxSearch(P_elem'*residual);
    % Bring the points in the order the've been chosen
    pts_elem = pts_elem(pts_local);
    % Add to global point set
    pts = [pts pts_elem];%#ok
    % Add points to basis (index is residual-dependent, so transfer to
    % global index)
    basis_idx = [basis_idx resi_idx(basis_local)];%#ok
    % Remove selected new basis from residual index set
    resi_idx(basis_local) = [];
    
    if isempty(resi_idx)
        aborted = true;
        break;
    end
    
    % Update residual
    basis = u(:,basis_idx);
    npts = length(pts);
    P = sparse(pts,1:npts,ones(npts,1),n,npts);
    c = (P'*basis) \ (P'*u(:,resi_idx));
    residual = u(:,resi_idx) - basis*c;
	
%     [a,b,c] = svd(Psel'*hlp,'econ');
%     basis = [basis o.orthonormalize(u*c)];%#ok
%     pts = [pts newpts];%#ok
%     npts = length(pts);
%     P = sparse(pts,1:npts,ones(npts,1),n,npts);
%     c = (P'*basis) \ (P'*u);
%     residual = u - basis*c;
%     hlp = u - basis*basis'*u;
end
elems = elems(1:elem_idx);

% no = sqrt(sigma*(residual.^2));
% maxerr(num_elems+1) = max(max(no,[],2));
% maxerr_res(num_elems+1) = max(abs(residual(:)));

pts_FEDEIM = pts;
err_FEDEIM = fd.getInterpolErrors(Uu, pts_FEDEIM, dof_vp_unass);
err_UDEIM_to_DEIM = fd.getInterpolErrors(Uu, pts_UDEIM, dof_vp_unass, s.f.Sigma);
elems_UDEIM = zeros(1,length(pts_UDEIM));
for k = 1:length(pts_UDEIM)
    elems_UDEIM(k) = sum(sum(s.f.Sigma(:,pts_UDEIM(1:k))*idx_elems(:,pts_UDEIM(1:k))') ~= 0);
end