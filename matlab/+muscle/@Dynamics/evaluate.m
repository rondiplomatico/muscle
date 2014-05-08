function duvw = evaluate(this, uvwdof, t)
    sys = this.System;
    mc = sys.Model.Config;
    g = mc.Geometry;
    fe_pos = mc.PosFE;
    fe_press = mc.PressFE;
    globidx_disp = sys.globidx_displ;
    globidx_press = sys.globidx_pressure;

    dofs_displ = fe_pos.NumNodes*3;

    % Cache variables instead of accessing them via this. in loops
    b1 = this.b1;
    d1 = this.d1;
    lfopt = this.lambdafopt;
    Pmax = this.Pmax;
    alpha = this.alpha; %#ok<*PROP>
    havefibres = sys.HasFibres;

    % Include dirichlet values to state vector
    uvwcomplete = zeros(2*dofs_displ + fe_press.NumNodes,1);
    uvwcomplete(sys.dof_idx_global) = uvwdof;
    uvwcomplete(sys.bc_dir_idx) = sys.bc_dir_val;

    % Init duv
    duvw = zeros(size(uvwcomplete));
    % THIS ALREADY FULFILLS THE u' = v ODE PART!
    duvw(1:dofs_displ) = uvwcomplete(dofs_displ+1:2*dofs_displ);

    dofsperelem_displ = fe_pos.DofsPerElement;
    dofsperelem_press = fe_press.DofsPerElement;
    num_gausspoints = g.NumGaussp;
    num_elements = fe_pos.NumElems;
    for m = 1:num_elements
        elemidx_displ = globidx_disp(:,:,m);
        elemidx_pressure = globidx_press(:,m);

        integrand_displ = zeros(3,dofsperelem_displ);
        integrand_press = zeros(dofsperelem_press,1);
        for gp = 1:num_gausspoints

            % Evaluate the pressure at gauss points
            w = uvwcomplete(elemidx_pressure);
            p = w' * fe_press.Ngp(:,gp,m);

            pos = 3*(gp-1)+1:3*gp;
            dtn = fe_pos.transgrad(:,pos,m);

            % Get coefficients for nodes of current element
            u = uvwcomplete(elemidx_displ);
            % Deformation gradient
            F = u * dtn;

            %% Isotropic part (Invariant I1 related)
            I1 = sum(sum((u'*u) .* (dtn*dtn')));
            
            %% Compile tensor
            P = p*inv(F)' + 2*(this.c10 + I1*this.c01)*F ...
                - 2*this.c01*F*(F'*F);

            %% Anisotropic part (Invariant I4 related)
            if havefibres
                dtna0 = sys.dtna0(:,gp,m);
                lambdafsq = sum(sum((u'*u) .* (dtna0*dtna0')));
                lambdaf = sqrt(lambdafsq);

                % Evaluate g function
                % Using a subfunction is 20% slower!
                % So: direct implementation here
                ratio = lambdaf/lfopt;
                fl = (-6.25*ratio*ratio + 12.5*ratio - 5.25) * (ratio >= .6) * (ratio <= 1.4);
                gval = (b1/lambdafsq)*(lambdaf^d1-1) + (Pmax/lambdaf)*fl*alpha;
                a0 = sys.a0oa0(:,:,(m-1)*num_gausspoints + gp);
                P = P + gval*F*a0;
            end
            
            weight = g.gaussw(gp) * fe_pos.elem_detjac(m,gp);

            integrand_displ = integrand_displ + weight * P * dtn';

            integrand_press = integrand_press + weight * (det(F)-1) * fe_press.Ngp(:,gp,m);
        end
        % We have v' + K(u) = 0, so the values of K(u) must be
        % written at the according locations of v'; those are the
        % same but +fielddofs later each
        outidx = elemidx_displ + dofs_displ;
        % Have MINUS here as the equation satisfies Mu'' + K(u,w) =
        % 0, but KerMor implements Mu'' = -K(u,w)
        duvw(outidx) = duvw(outidx) - integrand_displ;

        % Update pressure value at according positions
        duvw(elemidx_pressure) = duvw(elemidx_pressure) + integrand_press;
    end
    % Remove values at dirichlet nodes
    duvw(sys.bc_dir_idx) = [];

    if sys.UseDirectMassInversion
        duvw(sys.dof_idx_velo) = sys.Minv * duvw(sys.dof_idx_velo);
    end
end