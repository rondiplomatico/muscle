function duvw = evaluate(this, uvwdof, t)
    sys = this.System;
    mc = sys.Model.Config;
    fe_pos = mc.PosFE;
    geo = fe_pos.Geometry;
    fe_press = mc.PressFE;
    pgeo = fe_press.Geometry;
    globidx_disp = sys.globidx_displ;
    globidx_press = sys.globidx_pressure;

    dofs_pos = geo.NumNodes*3;

    % Cache variables instead of accessing them via this. in loops
    b1 = this.b1;
    d1 = this.d1;
    lfopt = this.lambdafopt;
    Pmax = this.Pmax;
    flfun = this.ForceLengthFun;
    alphaconst = this.alpha(t);
    havefibres = sys.HasFibres;
    havefibretypes = havefibres && ~isempty(mc.Pool);
    
%     visc = this.mu(1);
    if havefibretypes
        fibretypeweights = mc.FibreTypeWeights;
        % Input data is x1: fibre type, x2: mean current, x3: time
%         forceargs = [this.muprep; t*ones(1,size(this.muprep,2))];
        % This is the learned 
%         FibreForces = this.APExp.evaluate(forceargs)';
        FibreForces = mc.Pool.getActivation(t);
%         FibreForces = alphaconst*ones(size(this.muprep,2),1);
    end

    % Include dirichlet values to state vector
    uvwcomplete = zeros(2*dofs_pos + pgeo.NumNodes,1);
    uvwcomplete(sys.dof_idx_global) = uvwdof;
    uvwcomplete(sys.bc_dir_idx) = sys.bc_dir_val;
    % Check if velocity bc's should still be applied
    if t > sys.ApplyVelocityBCUntil
        uvwcomplete(sys.bc_dir_velo_idx) = 0;
    end

    % Init duvw
    duvw = zeros(size(uvwcomplete));
    % THIS ALREADY FULFILLS THE u' = v ODE PART!
    duvw(1:dofs_pos) = uvwcomplete(dofs_pos+1:2*dofs_pos);

    dofsperelem_pos = geo.DofsPerElement;
    dofsperelem_press = pgeo.DofsPerElement;
    num_gausspoints = fe_pos.GaussPointsPerElem;
    num_elements = geo.NumElements;
    for m = 1:num_elements
        elemidx_pos = globidx_disp(:,:,m);
        elemidx_velo = elemidx_pos+dofs_pos;
        elemidx_pressure = globidx_press(:,m);
        
        if havefibretypes 
            ftwelem = fibretypeweights(:,:,m)*FibreForces;
        end

        integrand_pos = zeros(3,dofsperelem_pos);
        integrand_press = zeros(dofsperelem_press,1);
        for gp = 1:num_gausspoints

            % Evaluate the pressure at gauss points
            w = uvwcomplete(elemidx_pressure);
            p = w' * fe_press.Ngp(:,gp,m);

            pos = 3*(gp-1)+1:3*gp;
            dtn = fe_pos.transgrad(:,pos,m);

            % Get coefficients for nodes of current element
            u = uvwcomplete(elemidx_pos);
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
%                 dtna0 = sys.dtna0(:,gp,m);
%                 lambdafsq = sum(sum((u'*u) .* (dtna0*dtna0')));
                lambdafsq = sum((F*sys.a0(:,gp,m)).^2);
                lambdaf = sqrt(lambdafsq);

                % Evaluate g function
                % Using a subfunction is 20% slower!
                % So: direct implementation here
                fl = flfun(lambdaf/lfopt);
                alpha = alphaconst;
                if havefibretypes 
                    alpha = ftwelem(gp);
                end
                markert = 0;
                if lambdaf > 1
                    markert = (b1/lambdafsq)*(lambdaf^d1-1);
                end
                gval = markert + (Pmax/lambdaf)*fl*alpha;
                a0 = sys.a0oa0(:,:,(m-1)*num_gausspoints + gp);
                P = P + gval*F*a0;
            end
            
%             Viscosity
%             if visc > 0
%                 v = uvwcomplete(elemidx_velo);
%                 P = P + visc * v * dtn;
% %                 P = P + visc * 2 * F' * (v * dtn);
%             end
            
            weight = fe_pos.GaussWeights(gp) * fe_pos.elem_detjac(m,gp);

            integrand_pos = integrand_pos + weight * P * dtn';

            integrand_press = integrand_press + weight * (det(F)-1) * fe_press.Ngp(:,gp,m);
        end
        % We have v' + K(u) = 0, so the values of K(u) must be
        % written at the according locations of v', i.e. elemidx_velo
        %
        % Have MINUS here as the equation satisfies Mu'' + K(u,w) =
        % 0, but KerMor implements Mu'' = -K(u,w)
        duvw(elemidx_velo) = duvw(elemidx_velo) - integrand_pos;

        % Update pressure value at according positions
        duvw(elemidx_pressure) = duvw(elemidx_pressure) + integrand_press;
    end
    
    %% Save & remove values at dirichlet pos/velo nodes
    this.LastBCResiduals = duvw([sys.bc_dir_displ_idx+dofs_pos; sys.bc_dir_velo_idx]);
    duvw(sys.bc_dir_idx) = [];
    
    fo = sum(duvw);
    fprintf('t=%20g, force=%20g\n',t,fo);
%     if fo > 1e4
%         keyboard;
%     end

    %% If direct mass matrix inversion is used
    if this.usemassinv
        duvw(sys.dof_idx_velo) = sys.Minv * duvw(sys.dof_idx_velo);
    end
end