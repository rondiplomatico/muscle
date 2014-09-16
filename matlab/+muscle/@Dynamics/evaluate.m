function duvw  = evaluate(this, uvwdof, t)
    sys = this.System;
    mc = sys.Model.Config;
    fe_pos = mc.PosFE;
    geo = fe_pos.Geometry;
    fe_press = mc.PressFE;
    pgeo = fe_press.Geometry;
    elem_idx_u_glob = sys.idx_u_glob_elems;
    elem_idx_p_glob = sys.idx_p_glob_elems;
    unassembled = this.ComputeUnassembled;

    num_u_glob = geo.NumNodes*3;
    num_v_glob = num_u_glob;

    % Cache variables instead of accessing them via this. in loops
    b1 = this.b1;
    d1 = this.d1;
    lfopt = this.lambdafopt;
    Pmax = this.Pmax;
    flfun = this.ForceLengthFun;
    alphaconst = this.alpha(t);
    havefibres = sys.HasFibres;
    havefibretypes = havefibres && size(uvwdof,1) > sys.num_uvp_dof;
    usecrossfibres = havefibres && this.crossfibres;
    if usecrossfibres
        b1cf = this.b1cf;
        d1cf = this.d1cf;
    end
    
%     visc = this.mu(1);
    if havefibretypes
        fibretypeweights = mc.FibreTypeWeights;
        % Input data is x1: fibre type, x2: mean current, x3: time
%         forceargs = [this.muprep; t*ones(1,size(this.muprep,2))];
        % This is the learned 
%         FibreForces = this.APExp.evaluate(forceargs)';
%         FibreForces = mc.Pool.getActivation(t);
        FibreForces = uvwdof(sys.num_uvp_dof+1:end);
%         FibreForces = alphaconst*ones(size(this.muprep,2),1);
    end

    % Include dirichlet values to state vector
    uvwcomplete = zeros(2*num_u_glob + pgeo.NumNodes,1);
    uvwcomplete(sys.idx_uv_dof_glob) = uvwdof(1:sys.num_uvp_dof);
    uvwcomplete(sys.idx_uv_bc_glob) = sys.val_uv_bc_glob;
    % Check if velocity bc's should still be applied
    if t > sys.ApplyVelocityBCUntil
        uvwcomplete(sys.idx_v_bc_glob) = 0;
    end
    
    dofsperelem_u = geo.DofsPerElement;
    dofsperelem_p = pgeo.DofsPerElement;
    num_gp = fe_pos.GaussPointsPerElem;
    num_elements = geo.NumElements;

    % Init result vector duvw
    if unassembled
        % dofs_pos for u', elems*3*dofperelem for v',
        % elems*dofsperelem_press for p'
        duvw = zeros(this.num_uvp_dof_unass,1);
        % Set u' = v
        duvw(1:num_u_glob) = uvwcomplete(num_u_glob + (1:num_v_glob));
        unass_offset_dvelo = num_u_glob;
        unass_offset_dpressure = unass_offset_dvelo + num_elements*3*dofsperelem_u;
    else
        duvw = zeros(size(uvwcomplete));
        % Set u' = v
        duvw(1:num_u_glob) = uvwcomplete(num_u_glob + (1:num_v_glob));
    end
    
    for m = 1:num_elements
        elemidx_u = elem_idx_u_glob(:,:,m); % 1:num_u_glob is all u
        elemidx_v = num_u_glob + elemidx_u; % num_u_glob next ones are all v
        elemidx_p = elem_idx_p_glob(:,m);
        
        u = uvwcomplete(elemidx_u);
        w = uvwcomplete(elemidx_p);
        
        if havefibretypes 
            ftwelem = fibretypeweights(:,:,m)*FibreForces;
        end

        integrand_u = zeros(3,dofsperelem_u);
        integrand_p = zeros(dofsperelem_p,1);
        for gp = 1:num_gp

            % Evaluate the pressure at gauss points
            p = w' * fe_press.Ngp(:,gp,m);

            pos = 3*(gp-1)+1:3*gp;
            dtn = fe_pos.transgrad(:,pos,m);

            % Deformation gradient
            F = u * dtn;
            C = F'*F;
           
            %% Isotropic part (Invariant I1 related)
%             I1 = sum(sum((u'*u) .* (dtn*dtn')));
            I1 = C(1,1) + C(2,2) + C(3,3);
            
            %% Compile tensor
            P = p*inv(F)' + 2*(this.c10 + I1*this.c01)*F ...
                - 2*this.c01*F*C;
            
            %% Anisotropic part (Invariant I4 related)
            if havefibres
                fibrenr = (m-1)*num_gp + gp;
                fibres = sys.a0Base(:,:,fibrenr);
                lambdaf = norm(F*fibres(:,1));
                
                % Evaluate g function
                % Using a subfunction is 20% slower!
                % So: direct implementation here
                fl = flfun(lambdaf/lfopt);
                alpha = alphaconst;
                if havefibretypes
                    alpha = ftwelem(gp);
                end
                markert = 0;
                % Using > 1 is deadly. All lambdas are equal to one at t=0
                % (reference config, analytical), but numerically this is
                % dependent on how precise F and hence lambda is computed.
                % It is very very close to one, but sometimes 1e-7 smaller
                % or bigger.. and that makes all the difference!
                if lambdaf > .999
                    markert = (b1/lambdaf^2)*(lambdaf^d1-1);
                end
                gval = markert + (Pmax/lambdaf)*fl*alpha; %+ 1000*max(0,(lambdaf-1))*alpha;
                P = P + gval*F*sys.a0oa0(:,:,fibrenr);
                
                %% Cross-fibre stiffness part
                if usecrossfibres
                    lambdaf = norm(F*fibres(:,2));
                    if lambdaf > .999
                        g1 = (b1cf/lambdaf^2)*(lambdaf^d1cf-1);
                        P = P + g1*F*sys.a0oa0n1(:,:,fibrenr);
                    end
                    lambdaf = norm(F*fibres(:,3));
                    if lambdaf > .999
                        g2 = (b1cf/lambdaf^2)*(lambdaf^d1cf-1);
                        P = P + g2*F*sys.a0oa0n2(:,:,fibrenr);
                    end
                end
            end
            
%             Viscosity
%             if visc > 0
%                 v = uvwcomplete(elemidx_velo);
%                 P = P + visc * v * dtn;
%             end

            weight = fe_pos.GaussWeights(gp) * fe_pos.elem_detjac(m,gp);

            integrand_u = integrand_u + weight * P * dtn';

            integrand_p = integrand_p + weight * (det(F)-1) * fe_press.Ngp(:,gp,m);
        end
        
        % Unassembled or assembled?
        if unassembled
            pos = unass_offset_dvelo + (1:3*dofsperelem_u) + (m-1) * 3 * dofsperelem_u;
            duvw(pos) = -integrand_u(:);
            pos = unass_offset_dpressure + (1:dofsperelem_p) + (m-1) * dofsperelem_p;
            duvw(pos) = integrand_p(:);
        else
            % We have v' + K(u) = 0, so the values of K(u) must be
            % written at the according locations of v', i.e. elemidx_velo
            %
            % Have MINUS here as the equation satisfies Mu'' + K(u,w) =
            % 0, but KerMor implements Mu'' = -K(u,w)
            duvw(elemidx_v) = duvw(elemidx_v) - integrand_u;
            duvw(elemidx_p) = duvw(elemidx_p) + integrand_p;
        end
    end
    
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