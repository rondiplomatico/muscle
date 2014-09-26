function [J, Jalpha, JLamDot] = getStateJacobian(this, uvwdof, t)
%     J = this.getStateJacobianFD(uvwdof,t);
%     return;
    sys = this.System;
    mc = sys.Model.Config; 
    fe_pos = mc.PosFE;
    geo = fe_pos.Geometry;
    fe_press = mc.PressFE;
    pgeo = fe_press.Geometry;

    N = geo.NumNodes;
    M = pgeo.NumNodes;
    dofs_pos = 3*N;
    dofsperelem_pos = geo.DofsPerElement;
    dofsperelem_press = pgeo.DofsPerElement;
    num_gausspoints = fe_pos.GaussPointsPerElem;
    num_elements = geo.NumElements;
    
    % Cache variables instead of accessing them via "this." in loops
    b1 = this.b1;
    d1 = this.d1;
    lfopt = this.lambdafopt;
    Pmax = this.Pmax;

    havefibres = sys.HasFibres;
    havefibretypes = sys.HasFibreTypes;
    usecrossfibres = this.crossfibres;
    hasforceargument = sys.HasForceArgument;
    if usecrossfibres
        b1cf = this.b1cf;
        d1cf = this.d1cf;
    end
    ldotpos = this.lambda_dot_pos;
    haveldotpos = ~isempty(ldotpos);
    if haveldotpos
       ildot = [];
       jldot = [];
       sldot = [];
    end
    
    %% Precompute the size of i,j,s for speed
    numXYZDofs_pos = 3*dofsperelem_pos;
    relidx_pos = 1:numXYZDofs_pos;
    relidx_press = 1:dofsperelem_press;
    % 3x for grad_u K(u,v,w), 1x for Grad_w K(u,v,w), 3 x pressure
    numparts = (3+1)*numXYZDofs_pos + dofsperelem_press;
%     if visc > 0
%         numparts = numparts+3*numXYZDofs_pos;
%     end
    totalsize = dofs_pos + num_elements*num_gausspoints*numparts;
    i = zeros(totalsize,1);
    j = zeros(totalsize,1);
    s = zeros(totalsize,1);
    
    %% Extra stuff if fibres are used
    if havefibretypes
        alphaconst = [];
        fibretypeweights = mc.FibreTypeWeights;
        nfibres = size(fibretypeweights,2);
        if sys.HasMotoPool
            FibreForces = mc.Pool.getActivation(t);
        elseif hasforceargument
            FibreForces = uvwdof(sys.num_uvp_dof+1:end) * min(1,t);
            totalsizeS = num_elements*num_gausspoints*nfibres;
            iS = zeros(totalsizeS,1);
            jS = zeros(totalsizeS,1);
            sS = zeros(totalsizeS,1);
            cur_offS = 0;
            columns_sarco_link = 53:56:56*nfibres;
        else
            error('No implemented');
        end
%         FibreForces = this.APExp.evaluate(forceargs)';
    else
        alphaconst = this.alpha(t);
    end

    %% -I part in u'(t) = -v(t)
    i(1:dofs_pos) = (1:dofs_pos)';
    j(1:dofs_pos) = ((1:dofs_pos)+dofs_pos)';
    s(1:dofs_pos) = 1;
    cur_off = dofs_pos;
    
    globidx_pos = sys.idx_u_glob_elems;
    globidx_press = sys.idx_p_glob_elems;

    % Include dirichlet values to state vector
    uvwcomplete = zeros(2*dofs_pos + pgeo.NumNodes,1);
    uvwcomplete(sys.idx_uv_dof_glob) = uvwdof(1:sys.num_uvp_dof);
    uvwcomplete(sys.idx_uv_bc_glob) = sys.val_uv_bc_glob;
    
    for m = 1:num_elements
        elemidx_u = globidx_pos(:,:,m);
        elemidx_v = elemidx_u + dofs_pos;
        elemidx_velo_linear = elemidx_v(:);
        elemidx_velo_linear3 = [elemidx_velo_linear
                                elemidx_velo_linear
                                elemidx_velo_linear];
        elemidx_p = globidx_press(:,m);
        elemidx_pressure3 = [elemidx_p
                             elemidx_p
                             elemidx_p];
        
        if havefibretypes 
            ftwelem = fibretypeweights(:,:,m)*FibreForces;
        end

        for gp = 1:num_gausspoints
            pos = 3*(gp-1)+1:3*gp;
            dtn = fe_pos.transgrad(:,pos,m);
            u = uvwcomplete(elemidx_u);
            
            % Deformation gradient
            F = u * dtn;

            Finv = inv(F);
            C = F'*F;

            % Evaluate the pressure at gauss points
            w = uvwcomplete(elemidx_p);
            p = w' * fe_press.Ngp(:,gp,m);

            % Invariant I1
            I1 = C(1,1) + C(2,2) + C(3,3);

            if havefibres
                %% Anisotropic part
                fibrenr = (m-1)*num_gausspoints + gp;
                fibres = sys.a0Base(:,:,fibrenr);
                Fa0 = F*fibres(:,1);
                lambdaf = norm(Fa0);

                ratio = lambdaf/lfopt;
                fl = this.ForceLengthFun(ratio);
                dfl = this.ForceLengthFunDeriv(ratio);
                alpha = alphaconst;
                if havefibretypes
                    alpha = ftwelem(gp);
                end
                alpha_prefactor = (Pmax/lambdaf)*fl;
                g_value = alpha_prefactor*alpha;
                dg_dlam = (Pmax/lambdaf^2)*alpha*(dfl - fl);
                % Using > 1 is deadly. All lambdas are equal to one at t=0
                % (reference config, analytical), but numerically this is
                % dependent on how precise F and hence lambda is computed.
                % It is very very close to one, but sometimes 1e-7 smaller
                % or bigger.. and that makes all the difference!
                if lambdaf > .999
                    g_value = g_value + (b1/lambdaf^2)*(lambdaf^d1-1);
                    dg_dlam = dg_dlam + (b1/lambdaf^3)*((d1-2)*lambdaf^d1 + 2);
                end
                a0 = sys.a0oa0(:,:,fibrenr);
                
                %% Cross-fibre stiffness part
                if usecrossfibres
                    Fcrossf1 = F*fibres(:,2);
                    lambdaa0_n1 = norm(Fcrossf1);
                    if lambdaa0_n1 > .999
                        xfibre1 = (b1cf/lambdaa0_n1^2)*(lambdaa0_n1^d1cf-1);
                        dxfibre1_dlam = (b1cf/lambdaa0_n1^3)*((d1cf-2)*lambdaa0_n1^d1cf + 2);
                    end
                    Fcrossf2 = F*fibres(:,3);
                    lambdaa0_n2 = norm(Fcrossf2);
                    if lambdaa0_n2 > .999
                        xfibre2 = (b1cf/lambdaa0_n2^2)*(lambdaa0_n2^d1cf-1);
                        dxfibre2_dlam = (b1cf/lambdaa0_n2^3)*((d1cf-2)*lambdaa0_n2^d1cf + 2);
                    end
                    a0oa0_2 = sys.a0oa0n1(:,:,fibrenr);
                    a0oa0_3 = sys.a0oa0n2(:,:,fibrenr);
                end
                
                %% Check if change rate of lambda at a certain point should be tracked
                if haveldotpos
                    k = find(ldotpos(1,:) == m & ldotpos(2,:) == gp);
                    if ~isempty(k)
                        Fdot = uvwcomplete(elemidx_v) * dtn;
                        Fdota0 = Fdot*fibres(:,1);
                        % Also update the current lambda_dot so that
                        % further calls within the getStateJacobian of the
                        % fullmuscle.Dynamics have the correct values!
                        this.lambda_dot(k) = Fa0'*Fdota0/lambdaf;

                        %% Assemble dLdot / du[i] and dLdot / dv[i]
                        JLamDot = zeros(2*numXYZDofs_pos,1);
                        for eldof = 1:dofsperelem_pos
                            % U_i^k = e_i dyad dPhik from script
                            U1k = [dtn(eldof,:); 0 0 0; 0 0 0];
                            U2k = [0 0 0; dtn(eldof,:); 0 0 0];
                            U3k = [0 0 0; 0 0 0; dtn(eldof,:)];
                            idx = (eldof-1)*3 + (1:3);
                            idxv = idx + 3*dofsperelem_pos;
                            %% x
                            Ua0 = U1k*fibres(:,1);
                            % dLdot / du[i]_x
                            JLamDot(idx(1)) = JLamDot(idx(1)) + Ua0'*Fdota0;
                            % dLdot / dv[i]_x
                            JLamDot(idxv(1)) = JLamDot(idxv(1)) + Fa0'*Ua0;

                            %% y
                            Ua0 = U2k*fibres(:,1);
                            % dLdot / du[i]_x
                            JLamDot(idx(2)) = JLamDot(idx(2)) + Ua0'*Fdota0;
                            % dLdot / dv[i]_x
                            JLamDot(idxv(2)) = JLamDot(idxv(2)) + Fa0'*Ua0;

                            %% z
                            Ua0 = U3k*fibres(:,1);
                            % dLdot / du[i]_x
                            JLamDot(idx(3)) = JLamDot(idx(3)) + Ua0'*Fdota0;
                            % dLdot / dv[i]_x
                            JLamDot(idxv(3)) = JLamDot(idxv(3)) + Fa0'*Ua0;
                        end
                        ildot = [ildot; k*ones(2*numXYZDofs_pos,1)]; %#ok
                        jldot = [jldot; elemidx_u(:); elemidx_v(:)];%#ok
                        sldot = [sldot; JLamDot];%#ok
                    end
                end
            end

            %% Main node loop
            % Precompute weight
            weight = fe_pos.GaussWeights(gp) * fe_pos.elem_detjac(m, gp);

            for k = 1:dofsperelem_pos
                % U_i^k = e_i dyad dPhik from script
                U1k = [dtn(k,:); 0 0 0; 0 0 0];
                U2k = [0 0 0; dtn(k,:); 0 0 0];
                U3k = [0 0 0; 0 0 0; dtn(k,:)];

                dI1duk = 2*sum([1; 1; 1] * (dtn(k,:) * dtn') .* u, 2);
                fac1 = 2*dI1duk*this.c01;
                fac2 = 2*(this.c10 + I1*this.c01);

                %% grad_u K(u,v,w)
                % Recall: gradients from nabla K_{u,w} are
                % negative, as KerMor implements Mu'' = -K(u,v,w)
                % instead of Mu'' + K(u,v,w) = 0
                
                % Assign i index as whole for x,y,z (speed)
                i(cur_off + (1:3*numXYZDofs_pos)) = elemidx_velo_linear3;
                
                % xdim
                dFtF1 = U1k'*F + F'*U1k;
                dPx = -p * (Finv * U1k * Finv)'...
                    + fac1(1)*F + fac2*U1k...
                    -2*this.c01 * (U1k * C + F*dFtF1);%#ok<*MINV>
                if havefibres
                    dlambda_dux = Fa0'*U1k*fibres(:,1)/lambdaf;
                    dPx = dPx + (dg_dlam*dlambda_dux*F + g_value*U1k)*a0;
                    if usecrossfibres 
                        if lambdaa0_n1 > .999
                            dlambda_dux = Fcrossf1'*U1k*fibres(:,2)/lambdaa0_n1;
                            dPx = dPx + (dxfibre1_dlam*dlambda_dux*F + xfibre1*U1k)*a0oa0_2;
                        end
                        if lambdaa0_n2 > .999
                            dlambda_dux = Fcrossf2'*U1k*fibres(:,3)/lambdaa0_n2;
                            dPx = dPx + (dxfibre2_dlam*dlambda_dux*F + xfibre2*U1k)*a0oa0_3;
                        end
                    end
                end
                j(cur_off + relidx_pos) = elemidx_u(1,k);
                snew = -weight * dPx * dtn';
                s(cur_off + relidx_pos) = snew(:);
                cur_off = cur_off + numXYZDofs_pos;

                % ydim
                dFtF2 = U2k'*F + F'*U2k;
                dPy = -p * (Finv * U2k * Finv)'...
                    +fac1(2)*F + fac2*U2k...
                    -2*this.c01 * (U2k * C + F*dFtF2);
                if havefibres
                    dlambda_duy = Fa0'*U2k*fibres(:,1)/lambdaf;
                    dPy = dPy + (dg_dlam*dlambda_duy*F + g_value*U2k)*a0;
                    if usecrossfibres 
                        if lambdaa0_n1 > .999
                            dlambda_duy = Fcrossf1'*U2k*fibres(:,2)/lambdaa0_n1;
                            dPy = dPy + (dxfibre1_dlam*dlambda_duy*F + xfibre1*U2k)*a0oa0_2;
                        end
                        if lambdaa0_n2 > .999
                            dlambda_duy = Fcrossf2'*U2k*fibres(:,3)/lambdaa0_n2;
                            dPy = dPy + (dxfibre2_dlam*dlambda_duy*F + xfibre2*U2k)*a0oa0_3;
                        end
                    end
                end
                j(cur_off + relidx_pos) = elemidx_u(2,k);
                snew = -weight * dPy * dtn';
                s(cur_off + relidx_pos) = snew(:);
                cur_off = cur_off + numXYZDofs_pos;

                % zdim
                dFtF3 = U3k'*F + F'*U3k;
                dPz = -p * (Finv * U3k * Finv)'...
                    +fac1(3)*F + fac2*U3k ...
                    - 2*this.c01 * (U3k * C + F*dFtF3);
                if havefibres
                    dlambda_duz = Fa0'*U3k*fibres(:,1)/lambdaf;
                    dPz = dPz + (dg_dlam*dlambda_duz*F + g_value*U3k)*a0;
                    if usecrossfibres 
                        if lambdaa0_n1 > .999
                            dlambda_duz = Fcrossf1'*U3k*fibres(:,2)/lambdaa0_n1;
                            dPz = dPz + (dxfibre1_dlam*dlambda_duz*F + xfibre1*U3k)*a0oa0_2;
                        end
                        if lambdaa0_n2 > .999
                            dlambda_duz = Fcrossf2'*U3k*fibres(:,3)/lambdaa0_n2;
                            dPz = dPz + (dxfibre2_dlam*dlambda_duz*F + xfibre2*U3k)*a0oa0_3;
                        end
                    end
                end
                j(cur_off + relidx_pos) = elemidx_u(3,k);
                snew = -weight * dPz * dtn';
                s(cur_off + relidx_pos) = snew(:);
                cur_off = cur_off + numXYZDofs_pos;
                
                %% grad_v K(u,v,w)
                % Viscosity part
%                 if visc > 0
%                     i(cur_off + relidx_pos) = elemidx_velo_linear;
%                     j(cur_off + relidx_pos) = elemidx_velo_XYZ(1,k);
%                     snew = -weight * visc * e1_dyad_dPhik * dtn';
%                     s(cur_off + relidx_pos) = snew(:);
%                     cur_off = cur_off + numXYZDofs_pos;
%                     
%                     % ydim
%                     i(cur_off + relidx_pos) = elemidx_velo_linear;
%                     j(cur_off + relidx_pos) = elemidx_velo_XYZ(2,k);
%                     snew = -weight * visc * e2_dyad_dPhik * dtn';
%                     s(cur_off + relidx_pos) = snew(:);
%                     cur_off = cur_off + numXYZDofs_pos;
%                     
%                     % zdim
%                     i(cur_off + relidx_pos) = elemidx_velo_linear;
%                     j(cur_off + relidx_pos) = elemidx_velo_XYZ(3,k);
%                     snew = -weight * visc * e3_dyad_dPhik * dtn';
%                     s(cur_off + relidx_pos) = snew(:);
%                     cur_off = cur_off + numXYZDofs_pos;
%                 end

                %% grad u g(u)
                precomp = weight * det(F) * fe_press.Ngp(:,gp,m);
                % Assign i index as whole for x,y,z (speed)
                i(cur_off + (1:3*dofsperelem_press)) = elemidx_pressure3;
                % dx
                j(cur_off + relidx_press) = elemidx_u(1,k);
                s(cur_off + relidx_press) = sum(diag(Finv*U1k)) * precomp;
                cur_off = cur_off + dofsperelem_press;
                % dy
                j(cur_off + relidx_press) = elemidx_u(2,k);
                s(cur_off + relidx_press) = sum(diag(Finv*U2k)) * precomp;
                cur_off = cur_off + dofsperelem_press;
                %dz
                j(cur_off + relidx_press) = elemidx_u(3,k);
                s(cur_off + relidx_press) = sum(diag(Finv*U3k)) * precomp;
                cur_off = cur_off + dofsperelem_press;
            end
            
            %% Grad_w K(u,w)
            for k = 1:dofsperelem_press
                i(cur_off + relidx_pos) = elemidx_velo_linear;
                j(cur_off + relidx_pos) = elemidx_p(k);
                snew = -weight * fe_press.Ngp(k,gp,m) * Finv' * dtn';
                s(cur_off + relidx_pos) = snew(:);
                cur_off = cur_off + numXYZDofs_pos;
            end
            
            %% Grad_s K(u,w,s)
            if hasforceargument
                for k = 1:nfibres
                    dPsk = alpha_prefactor * fibretypeweights(gp,k,m) * F * a0;
                    iS(cur_offS + relidx_pos) = elemidx_velo_linear-dofs_pos;
                    jS(cur_offS + relidx_pos) = columns_sarco_link(k);
                    snew = -weight * dPsk * dtn';
                    sS(cur_offS + relidx_pos) = snew(:);
                    cur_offS = cur_offS + numXYZDofs_pos;
                end
            end
        end
    end
    J = sparse(i,j,s,6*N+M,6*N+M);
    % Remove values at dirichlet nodes
    J(:,sys.idx_uv_bc_glob) = [];
    J(sys.idx_uv_bc_glob,:) = [];
    
    if this.usemassinv
        % Multiply with inverse of Mass matrix!
        J(sys.idx_v_dof_glob,:) = sys.Minv*J(sys.idx_v_dof_glob,:);
    end
    
    Jalpha = [];
    if hasforceargument
        Jalpha = sparse(iS,jS,sS,3*N,nfibres*56);
        % Remove those that are connected to dirichlet values
        Jalpha([sys.idx_u_bc_glob; sys.idx_v_bc_glob],:) = [];
    end
    
    JLamDot = [];
    if haveldotpos
        JLamDot = sparse(ildot,double(jldot),sldot,this.nfibres,6*N);
        JLamDot(:,sys.idx_uv_bc_glob) = [];
    end
end