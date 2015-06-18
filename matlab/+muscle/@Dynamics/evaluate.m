function duvw  = evaluate(this, uvwdof, t)
    % This function represents the entire nonlinear part of the transformed
    % first order system.
    %
    % It performs the substitution u'=v and invokes the nonlinear stiffness
    % operator K and the algebraic constraint operator g(u)
    %
    
    this.nfevals = this.nfevals+1;
    
    sys = this.System;
    m = sys.Model;
    mc = m.Config;
    fe_pos = mc.PosFE;
    geo = fe_pos.Geometry;
    fe_press = mc.PressFE;
    pgeo = fe_press.Geometry;
    unassembled = this.ComputeUnassembled;

    num_u_glob = geo.NumNodes*3;
    num_v_glob = num_u_glob;
    isproj = ~isempty(this.V);
    
    % If we evaluate inside a projected (reduced) model, reconstruct 
    if isproj
        % If projection is done before transformation to the first order
        % system, we can simply forward the coefficients for reduced v to
        % the coefficients for reduced u!
        if m.ProjectionFirst
            effsize = m.Data.ProjectionSpaces(1).LastEffectiveSize;
            
            duvw = zeros(size(uvwdof));
            duvw(1:effsize) = uvwdof((effsize+1):2*effsize);
        end
        
        % In any case: reconstruct full uvw dof vector, as the K operator
        % might need u,v,w quantities
        uvwdof = this.V*uvwdof;
    end

    %% Include dirichlet values to state vector
    uvwcomplete = zeros(2*num_u_glob + pgeo.NumNodes,1);
    % Insert DOFs
    uvwcomplete(sys.idx_uv_dof_glob) = uvwdof(1:sys.num_uvp_dof);
    % Insert BC
    uvwcomplete(sys.idx_uv_bc_glob) = sys.val_uv_bc_glob;
    % Check if velocity bc's should be applied in time-dependent manner
    if ~isempty(this.velo_bc_fun)
        uvwcomplete(sys.idx_v_bc_glob) = ...
            this.velo_bc_fun(t)*uvwcomplete(sys.idx_v_bc_glob);
    end
    
    if ~(isproj && m.ProjectionFirst)
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
    end
    
    %% Evaluate K(u,v,w) and g(u)
    % This is the main FEM-loop
    dvw = this.Kg(uvwcomplete,t);
    
    if ~(isproj && m.ProjectionFirst)
        duvw((num_u_glob+1):end) = dvw;
        
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
    end

    if isproj
        if m.ProjectionFirst
            dvw(sys.idx_u_bc_glob) = [];
            Wpart = this.W((sys.num_u_dof+1):end,(effsize+1):end);
            duvw((effsize+1):end) = Wpart'*dvw;
        else
            duvw = this.W'*duvw;
        end
    end
    
end