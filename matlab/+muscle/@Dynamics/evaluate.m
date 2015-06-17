function duvw  = evaluate(this, uvwdof, t)
    % This function represents the entire nonlinear part of the transformed
    % first order system.
    %
    % It performs the substitution u'=v and invokes the nonlinear stiffness
    % operator K and the algebraic constraint operator g(u)
    %
    
    this.nfevals = this.nfevals+1;
    
    sys = this.System;
    mc = sys.Model.Config;
    fe_pos = mc.PosFE;
    geo = fe_pos.Geometry;
    fe_press = mc.PressFE;
    pgeo = fe_press.Geometry;
    unassembled = this.ComputeUnassembled;

    num_u_glob = geo.NumNodes*3;
    num_v_glob = num_u_glob;
    
    % If we evaluate inside a projected (reduced) model, reconstruct 
    if ~isempty(this.V)
        uvwdof = this.V*uvwdof;
    end

    %% Include dirichlet values to state vector
    uvwcomplete = zeros(2*num_u_glob + pgeo.NumNodes,1);
    uvwcomplete(sys.idx_uv_dof_glob) = uvwdof(1:sys.num_uvp_dof);
    uvwcomplete(sys.idx_uv_bc_glob) = sys.val_uv_bc_glob;
    % Check if velocity bc's should be applied in time-dependent manner
    if ~isempty(this.velo_bc_fun)
        uvwcomplete(sys.idx_v_bc_glob) = ...
            this.velo_bc_fun(t)*uvwcomplete(sys.idx_v_bc_glob);
    end
    
    % Init result vector duvw
    if unassembled
        % dofs_pos for u', elems*3*dofperelem for v',
        % elems*dofsperelem_press for p'
        duvw = zeros(this.num_uvp_dof_unass,1);    
    else
        duvw = zeros(size(uvwcomplete));
    end
    
    %% Set u' = v
    duvw(1:num_u_glob) = uvwcomplete(num_u_glob + (1:num_v_glob));
    
    %% Evaluate K(u,v,w) and g(u)
    % This is the main FEM-loop
    duvw((num_u_glob+1):end) = this.Kg(uvwcomplete,t);
    
    %% Remove dirichlet boundary condition values
    if unassembled
        duvw(this.idx_uv_bc_glob_unass) = [];
    else
        %% Save & remove values at dirichlet pos/velo nodes
        this.LastBCResiduals = duvw([sys.idx_u_bc_glob+num_u_glob; sys.idx_v_bc_glob]);
        duvw(sys.idx_uv_bc_glob) = [];
        
        %% If direct mass matrix inversion is used
        if this.usemassinv
            duvw(sys.idx_v_dof_glob) = sys.Minv * duvw(sys.idx_v_dof_glob);
        end
    end
    
    if ~isempty(this.W)
        duvw = this.W'*duvw;
    end
    
end