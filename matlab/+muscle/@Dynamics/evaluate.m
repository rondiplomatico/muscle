function dK  = evaluate(this, uvwdof, t)
    % This function represents the entire nonlinear part of the transformed
    % first order system.
    %
    % It performs the substitution u'=v and invokes the nonlinear stiffness
    % operator K and the algebraic constraint operator g(u)
    %
    
    this.nfevals = this.nfevals+1;
    
    sys = this.System;
%     m = sys.Model;
%     mc = m.Config;
%     fe_pos = mc.PosFE;
%     geo = fe_pos.Geometry;
%     fe_press = mc.PressFE;
%     pgeo = fe_press.Geometry;
%     unassembled = this.ComputeUnassembled;

%     num_u_glob = geo.NumNodes*3;
%     num_v_glob = num_u_glob;
%     isproj = ~isempty(this.V);
    
%     % If we evaluate inside a projected (reduced) model, reconstruct 
%     if isproj
%         effsize_reduced_u_dofs = this.reduced_space_size;
% 
%         % Set z'=w directly
%         duvw = zeros(size(uvwdof));
%         hlp = 2*effsize_reduced_u_dofs;
%         duvw(1:effsize_reduced_u_dofs) = ...
%             uvwdof((effsize_reduced_u_dofs+1):hlp);
%         
%         % Include velocity boundary conditions - they are not-projected
%         % dofs coming just after the actual effsize_reduced_u_dofs.
%         if sys.num_v_bc > 0
%             velo_bc = sys.val_v_bc;
%             if ~isempty(this.velo_bc_fun)
%                 velo_bc = this.velo_bc_fun(t)*velo_bc;
%             end
%             duvw(hlp+(1:sys.num_v_bc)) = velo_bc;
%         end
%         
%         % Then: reconstruct full uvw dof vector, as the K operator
%         % might need u,v,w quantities
%         uvwdof = this.V*uvwdof;
%     end

    %% Include dirichlet values to state vector
    uvwcomplete = sys.includeDirichletValues(t, uvwdof);
    
%     uvwcomplete = zeros(2*num_u_glob + pgeo.NumNodes,1);
%     % Insert DOFs
%     uvwcomplete(sys.idx_uv_dof_glob) = uvwdof(1:sys.NumTotalDofs);
%     % Insert BC
%     uvwcomplete(sys.idx_uv_bc_glob) = sys.val_uv_bc_glob;
    
%     if ~isproj
%         % Init result vector duvw
%         if unassembled
%             % dofs_pos for u', elems*3*dofperelem for v',
%             % elems*dofsperelem_press for p'
%             duvw = zeros(this.NumTotalDofs_unass,1);    
%         else
%             dK = zeros(size(uvwcomplete));
%         end
%     end
    
    %% Evaluate K(u,v,w) and g(u)
    % This is the main FEM-loop, which evaluates K(u,v,w) and g(u)
    % simultaneously (efficiency for only one FEM-loop is required)
    dKg = this.Kg(uvwcomplete,t);
    
    % Extract boundary condition residuals for later use
    this.LastBCResiduals = dKg(sys.idx_v_bc_local);
    dKg(sys.idx_v_bc_local) = [];
    
    dd = sys.NumDerivativeDofs;
    % Cache the algebraic constraint evaluation - will be read by the
    % ConstraintsFun. This is just an efficiency hack that avoids having to
    % run a second FEM loop, as the constraint condition can be computed
    % along.
    this.curGC = dKg(dd+1:end);
    dK = dKg(1:dd);
    
%     if isproj
%         % Kick out dirichlet values from full state space vector of Kg part
%         bc_idx = sys.idx_u_bc_glob;
%         if ~isempty(sys.idx_v_bc_glob)
%             bc_idx = [bc_idx; sys.idx_v_bc_glob-num_u_glob];
%         end
%         dvw(bc_idx) = [];
%         
%         % Compute the position of the Kg-dof part within the current
%         % reduced space
%         red_pos = (effsize_reduced_u_dofs+1):size(this.W,2);
%         if sys.num_v_bc > 0
%             red_pos(effsize_reduced_u_dofs+(1:sys.num_v_bc)) = [];
%         end
%         dKg(red_pos) = this.W((sys.NumStateDofs+1):end,red_pos)'*dvw;
%     else
%         dKg((num_u_glob+1):end) = dvw;
%         
%         %% Remove dirichlet boundary condition values
%         if unassembled
%             dKg(this.idx_uv_bc_glob_unass) = [];
%         else
%             %% Save & remove values at dirichlet pos/velo nodes
%             this.LastBCResiduals = dKg([sys.idx_u_bc_glob+num_u_glob; sys.idx_v_bc_glob]);
%             dKg(sys.idx_uv_bc_glob) = [];
% 
%             %% If direct mass matrix inversion is used
%             if this.usemassinv
%                 dKg(sys.idx_v_dof_glob) = sys.Minv * dKg(sys.idx_v_dof_glob);
%             end
%         end
%     end
end