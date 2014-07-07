function J = getStateJacobian(this, uvwdof, t)
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
%     visc = this.mu(1);
    alphaconst = this.alpha(t);
    havefibres = sys.HasFibres;
    havefibretypes = havefibres && ~isempty(mc.Pool);
    usecrossfibres = havefibres && this.crossfibres;
    if usecrossfibres
        b1cf = this.b1cf;
        d1cf = this.d1cf;
    end
    
    if havefibretypes
        fibretypeweights = mc.FibreTypeWeights;
        % Input data is x1: fibre type, x2: mean current, x3: time
%         forceargs = [this.muprep; t*ones(1,size(this.muprep,2))];
        % This is the learned 
%         FibreForces = this.APExp.evaluate(forceargs)';
%         FibreForces = alphaconst*ones(size(this.muprep,2),1);
        FibreForces = mc.Pool.getActivation(t);
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

    %% -I part in u'(t) = -v(t)
    i(1:dofs_pos) = (1:dofs_pos)';
    j(1:dofs_pos) = ((1:dofs_pos)+dofs_pos)';
    s(1:dofs_pos) = 1;
    cur_off = dofs_pos;
    
    globidx_pos = sys.globidx_displ;
    globidx_press = sys.globidx_pressure;

    % Include dirichlet values to state vector
    uvwcomplete = zeros(2*dofs_pos + pgeo.NumNodes,1);
    uvwcomplete(sys.dof_idx_global) = uvwdof;
    uvwcomplete(sys.bc_dir_idx) = sys.bc_dir_val;
    
    for m = 1:num_elements
        elemidx_pos_XYZ = globidx_pos(:,:,m);
        elemidx_velo_XYZ = elemidx_pos_XYZ + dofs_pos;
        elemidx_velo_linear = elemidx_velo_XYZ(:);
        elemidx_velo_linear3 = [elemidx_velo_linear
                                elemidx_velo_linear
                                elemidx_velo_linear];
        elemidx_pressure = globidx_press(:,m);
        elemidx_pressure3 = [elemidx_pressure
                             elemidx_pressure
                             elemidx_pressure];
        
        if havefibretypes 
            ftwelem = fibretypeweights(:,:,m)*FibreForces;
        end

        for gp = 1:num_gausspoints
            pos = 3*(gp-1)+1:3*gp;
            dtn = fe_pos.transgrad(:,pos,m);
            u = uvwcomplete(elemidx_pos_XYZ);
            
            % Deformation gradient
            F = u * dtn;

            Finv = inv(F);
            C = F'*F;

            % Evaluate the pressure at gauss points
            w = uvwcomplete(elemidx_pressure);
            p = w' * fe_press.Ngp(:,gp,m);

            % Invariant I1
            I1 = C(1,1) + C(2,2) + C(3,3);

            if havefibres
                %% Anisotropic part
                fibrenr = (m-1)*num_gausspoints + gp;
                fibres = sys.a0Base(:,:,fibrenr);
                Fa0 = F*fibres(:,1);
                lambda_a0 = norm(Fa0);

                ratio = lambda_a0/lfopt;
                fl = this.ForceLengthFun(ratio);
                dfl = this.ForceLengthFunDeriv(ratio);
                alpha = alphaconst;
                if havefibretypes
                    alpha = ftwelem(gp);
                end
                g_value = (Pmax/lambda_a0)*fl*alpha;
                dg_dlam = (Pmax/lambda_a0^2)*alpha*(dfl - fl);
                % Using > 1 is deadly. All lambdas are equal to one at t=0
                % (reference config, analytical), but numerically this is
                % dependent on how precise F and hence lambda is computed.
                % It is very very close to one, but sometimes 1e-7 smaller
                % or bigger.. and that makes all the difference!
                if lambda_a0 > .999
                    g_value = g_value + (b1/lambda_a0^2)*(lambda_a0^d1-1);
                    dg_dlam = dg_dlam + (b1/lambda_a0^3)*((d1-2)*lambda_a0^d1 + 2);
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
                    dlambda_dux = Fa0'*U1k*fibres(:,1)/lambda_a0;
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
                j(cur_off + relidx_pos) = elemidx_pos_XYZ(1,k);
                snew = -weight * dPx * dtn';
                s(cur_off + relidx_pos) = snew(:);
                cur_off = cur_off + numXYZDofs_pos;

                % ydim
                dFtF2 = U2k'*F + F'*U2k;
                dPy = -p * (Finv * U2k * Finv)'...
                    +fac1(2)*F + fac2*U2k...
                    -2*this.c01 * (U2k * C + F*dFtF2);
                if havefibres
                    dlambda_duy = Fa0'*U2k*fibres(:,1)/lambda_a0;
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
                j(cur_off + relidx_pos) = elemidx_pos_XYZ(2,k);
                snew = -weight * dPy * dtn';
                s(cur_off + relidx_pos) = snew(:);
                cur_off = cur_off + numXYZDofs_pos;

                % zdim
                dFtF3 = U3k'*F + F'*U3k;
                dPz = -p * (Finv * U3k * Finv)'...
                    +fac1(3)*F + fac2*U3k ...
                    - 2*this.c01 * (U3k * C + F*dFtF3);
                if havefibres
                    dlambda_duz = Fa0'*U3k*fibres(:,1)/lambda_a0;
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
                j(cur_off + relidx_pos) = elemidx_pos_XYZ(3,k);
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
                j(cur_off + relidx_press) = elemidx_pos_XYZ(1,k);
                s(cur_off + relidx_press) = sum(diag(Finv*U1k)) * precomp;
                cur_off = cur_off + dofsperelem_press;
                % dy
                j(cur_off + relidx_press) = elemidx_pos_XYZ(2,k);
                s(cur_off + relidx_press) = sum(diag(Finv*U2k)) * precomp;
                cur_off = cur_off + dofsperelem_press;
                %dz
                j(cur_off + relidx_press) = elemidx_pos_XYZ(3,k);
                s(cur_off + relidx_press) = sum(diag(Finv*U3k)) * precomp;
                cur_off = cur_off + dofsperelem_press;
            end
            %% Grad_w K(u,v,w)
            for k = 1:dofsperelem_press
                i(cur_off + relidx_pos) = elemidx_velo_linear;
                j(cur_off + relidx_pos) = elemidx_pressure(k);
                snew = -weight * fe_press.Ngp(k,gp,m) * Finv' * dtn';
                s(cur_off + relidx_pos) = snew(:);
                cur_off = cur_off + numXYZDofs_pos;
            end
        end
    end
    J = sparse(i,j,s,6*N+M,6*N+M);
    % Remove values at dirichlet nodes
    J(:,sys.bc_dir_idx) = [];
    J(sys.bc_dir_idx,:) = [];

    if this.usemassinv
        % Multiply with inverse of Mass matrix!
        J(sys.dof_idx_velo,:) = sys.Minv*J(sys.dof_idx_velo,:);
    end
end