function computeSparsityPattern(this)
    sys = this.System;
    mc = sys.Model.Config;
    fe_pos = mc.PosFE;
    geo = fe_pos.Geometry;
    fe_press = mc.PressFE;
    pgeo = fe_press.Geometry;

    N = geo.NumNodes;
    M = pgeo.NumNodes;

    %% -I part in u'(t) = -v(t)
    i = (1:3*N)';
    j = ((1:3*N)+3*N)';

    globidx_disp = sys.globidx_displ;
    globidx_press = sys.globidx_pressure;

    dofs_displ = N*3;
    visc = this.fViscosity;

    dofsperelem_displ = geo.DofsPerElement;
    dofsperelem_press = pgeo.DofsPerElement;
    num_gausspoints = geo.GaussPointsPerElem;
    num_elements = geo.NumElements;
    for m = 1:num_elements
        elemidx_displ = globidx_disp(:,:,m);
        elemidx_velo = elemidx_displ + dofs_displ;
        elemidx_pressure = globidx_press(:,m);
        inew = elemidx_velo(:);
        one = ones(size(inew));
        for gp = 1:num_gausspoints
            for k = 1:dofsperelem_displ
                %% Grad_u K(u,v,w)
                % xdim
                i = [i; inew]; %#ok<*AGROW>
                j = [j; one*elemidx_displ(1,k)];
                % ydim
                i = [i; inew]; 
                j = [j; one*elemidx_displ(2,k)]; 
                % zdim
                i = [i; inew]; 
                j = [j; one*elemidx_displ(3,k)]; 

                %% Grad_v K(u,v,w)
                if visc > 0
                    % xdim
                    i = [i; inew]; %#ok<*AGROW>
                    j = [j; one*elemidx_velo(1,k)];
                    % ydim
                    i = [i; inew]; 
                    j = [j; one*elemidx_velo(2,k)]; 
                    % zdim
                    i = [i; inew]; 
                    j = [j; one*elemidx_velo(3,k)]; 
                end

                %% grad u g(u)
                % dx
                i = [i; elemidx_pressure(:)];
                j = [j; ones(dofsperelem_press,1)*elemidx_displ(1,k)]; 
                % dy
                i = [i; elemidx_pressure(:)];
                j = [j; ones(dofsperelem_press,1)*elemidx_displ(2,k)]; 
                %dz
                i = [i; elemidx_pressure(:)];
                j = [j; ones(dofsperelem_press,1)*elemidx_displ(3,k)]; 
            end
            %% Grad_w K(u,v,w)
            inew = elemidx_velo(:);
            for k = 1:dofsperelem_press
                i = [i; inew];
                j = [j; ones(3*dofsperelem_displ,1)*elemidx_pressure(k)]; 
            end
        end
    end
    J = sparse(i,j,ones(size(i)),6*N+M,6*N+M);
    % Remove values at dirichlet nodes
    J(:,sys.bc_dir_idx) = [];
    J(sys.bc_dir_idx,:) = [];

    this.JSparsityPattern = logical(J);
    end

