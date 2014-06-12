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
        elemidx_pressure = globidx_press(:,m);
        
%         one = ones(numidx,1);
        
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
                a0pos = (m-1)*num_gausspoints + gp;
                a0B = sys.a0Base(:,:,a0pos);
                a0BI = sys.a0BaseInv(:,:,a0pos);
                lambdas = a0BI*F*a0B;
                lambda_a0 = lambdas(1,1);

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
                a0 = sys.a0oa0(:,:,(m-1)*num_gausspoints + gp);
                
                %% Cross-fibre stiffness part
                if usecrossfibres
                    lambda_a0n1 = lambdas(2,2);
                    if lambda_a0n1 > .999
                        xfibre1 = (b1cf/lambda_a0n1^2)*(lambda_a0n1^d1cf-1);
                        dxfibre1_dlam = (b1cf/lambda_a0n1^3)*((d1cf-2)*lambda_a0n1^d1cf + 2);
                    end
                    lambda_a0n2 = lambdas(3,3);
                    if lambda_a0n2 > .999
                        xfibre2 = (b1cf/lambda_a0n2^2)*(lambda_a0n2^d1cf-1);
                        dxfibre2_dlam = (b1cf/lambda_a0n2^3)*((d1cf-2)*lambda_a0n2^d1cf + 2);
                    end
                    
                    a0n1 = sys.a0oa0n1(:,:,(m-1)*num_gausspoints + gp);
                    a0n2 = sys.a0oa0n2(:,:,(m-1)*num_gausspoints + gp);
                end
            end

            %% Main node loop
            % Precompute weight
            weight = fe_pos.GaussWeights(gp) * fe_pos.elem_detjac(m, gp);

            for k = 1:dofsperelem_pos
                e1_dyad_dPhik = [dtn(k,:); 0 0 0; 0 0 0];
                e2_dyad_dPhik = [0 0 0; dtn(k,:); 0 0 0];
                e3_dyad_dPhik = [0 0 0; 0 0 0; dtn(k,:)];

                dI1duk = 2*sum([1; 1; 1] * (dtn(k,:) * dtn') .* u, 2);
                fac1 = 2*dI1duk*this.c01;
                fac2 = 2*(this.c10 + I1*this.c01);

                % correct! :-)
                if havefibres
                    dlamdas_x = diag(a0B'*e1_dyad_dPhik*a0B);
                    dlamdas_y = diag(a0B'*e2_dyad_dPhik*a0B);
                    dlamdas_z = diag(a0B'*e3_dyad_dPhik*a0B);
                end

                %% grad_u K(u,v,w)
                % Recall: gradients from nabla K_{u,w} are
                % negative, as KerMor implements Mu'' = -K(u,v,w)
                % instead of Mu'' + K(u,v,w) = 0
                % xdim
                dFtF1 = e1_dyad_dPhik'*F + F'*e1_dyad_dPhik;
                dPx = -p * (Finv * e1_dyad_dPhik * Finv)'...
                    + fac1(1)*F + fac2*e1_dyad_dPhik...
                    -2*this.c01 * (e1_dyad_dPhik * C + F*dFtF1);%#ok<*MINV>
                if havefibres
                    dPx = dPx + (dg_dlam*dlamdas_x(1)*F + g_value*e1_dyad_dPhik)*a0;
                    if usecrossfibres 
                        if lambda_a0n1 > .999
                            dPx = dPx + (dxfibre1_dlam*dlamdas_x(2)*F + xfibre1*e1_dyad_dPhik)*a0n1;
                        end
                        if lambda_a0n2 > .999
                            dPx = dPx + (dxfibre2_dlam*dlamdas_x(3)*F + xfibre2*e1_dyad_dPhik)*a0n2;
                        end
                    end
                end
                i(cur_off + relidx_pos) = elemidx_velo_linear;
                j(cur_off + relidx_pos) = elemidx_pos_XYZ(1,k);
                snew = -weight * dPx * dtn';
                s(cur_off + relidx_pos) = snew(:);
                cur_off = cur_off + numXYZDofs_pos;

                % ydim
                dFtF2 = e2_dyad_dPhik'*F + F'*e2_dyad_dPhik;
                dPy = -p * (Finv * e2_dyad_dPhik * Finv)'...
                    +fac1(2)*F + fac2*e2_dyad_dPhik...
                    -2*this.c01 * (e2_dyad_dPhik * C + F*dFtF2);
                if havefibres
                    dPy = dPy + (dg_dlam*dlamdas_y(1)*F + g_value*e2_dyad_dPhik)*a0;
                    if usecrossfibres 
                        if lambda_a0n1 > .999
                            dPy = dPy + (dxfibre1_dlam*dlamdas_y(2)*F + xfibre1*e2_dyad_dPhik)*a0n1;
                        end
                        if lambda_a0n2 > .999
                            dPy = dPy + (dxfibre2_dlam*dlamdas_y(3)*F + xfibre2*e2_dyad_dPhik)*a0n2;
                        end
                    end
                end
                i(cur_off + relidx_pos) = elemidx_velo_linear;
                j(cur_off + relidx_pos) = elemidx_pos_XYZ(2,k);
                snew = -weight * dPy * dtn';
                s(cur_off + relidx_pos) = snew(:);
                cur_off = cur_off + numXYZDofs_pos;

                % zdim
                dFtF3 = e3_dyad_dPhik'*F + F'*e3_dyad_dPhik;
                dPz = -p * (Finv * e3_dyad_dPhik * Finv)'...
                    +fac1(3)*F + fac2*e3_dyad_dPhik ...
                    - 2*this.c01 * (e3_dyad_dPhik * C + F*dFtF3);
                if havefibres
                    dPz = dPz + (dg_dlam*dlamdas_z(1)*F + g_value*e3_dyad_dPhik)*a0;
                    if usecrossfibres 
                        if lambda_a0n1 > .999
                            dPz = dPz + (dxfibre1_dlam*dlamdas_z(2)*F + xfibre1*e3_dyad_dPhik)*a0n1;
                        end
                        if lambda_a0n2 > .999
                            dPz = dPz + (dxfibre2_dlam*dlamdas_z(3)*F + xfibre2*e3_dyad_dPhik)*a0n2;
                        end
                    end
                end
                i(cur_off + relidx_pos) = elemidx_velo_linear;
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
                % dx
                i(cur_off + relidx_press) = elemidx_pressure;
                j(cur_off + relidx_press) = elemidx_pos_XYZ(1,k);
                s(cur_off + relidx_press) = sum(diag(Finv*e1_dyad_dPhik)) * precomp;
                cur_off = cur_off + dofsperelem_press;
%                 i = [i; elemidx_pressure];
%                 j = [j; ones(dofsperelem_press,1)*elemidx_pos_XYZ(1,k)];
%                 s = [s; sum(diag(Finv*e1_dyad_dPhik)) * precomp];
                % dy
                i(cur_off + relidx_press) = elemidx_pressure;
                j(cur_off + relidx_press) = elemidx_pos_XYZ(2,k);
                s(cur_off + relidx_press) = sum(diag(Finv*e2_dyad_dPhik)) * precomp;
                cur_off = cur_off + dofsperelem_press;
%                 i = [i; elemidx_pressure];
%                 j = [j; ones(dofsperelem_press,1)*elemidx_pos_XYZ(2,k)];
%                 s = [s; sum(diag(Finv*e2_dyad_dPhik)) * precomp];
                %dz
                i(cur_off + relidx_press) = elemidx_pressure;
                j(cur_off + relidx_press) = elemidx_pos_XYZ(3,k);
                s(cur_off + relidx_press) = sum(diag(Finv*e3_dyad_dPhik)) * precomp;
                cur_off = cur_off + dofsperelem_press;
%                 i = [i; elemidx_pressure];
%                 j = [j; ones(dofsperelem_press,1)*elemidx_pos_XYZ(3,k)];
%                 s = [s; sum(diag(Finv*e3_dyad_dPhik)) * precomp];
            end
            %% Grad_w K(u,v,w)
            for k = 1:dofsperelem_press
%                 i = [i; i_veloidx];
%                 j = [j; ones(3*dofsperelem_pos,1)*elemidx_pressure(k)];
%                 snew = -weight * fe_press.Ngp(k,gp,m) * Finv' * dtn';
%                 s = [s; snew(:)];
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