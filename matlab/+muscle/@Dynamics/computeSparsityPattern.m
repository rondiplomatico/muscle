function [SP, SPalpha, SPLamDot] = computeSparsityPattern(this)
    sys = this.System;
    mc = sys.Model.Config;
    fe_pos = mc.PosFE;
    geo = fe_pos.Geometry;
    fe_press = mc.PressFE;
    pgeo = fe_press.Geometry;

    N = geo.NumNodes;
    M = pgeo.NumNodes;
    
    % Jalpha
    if this.nfibres > 0
        iS = [];
        jS = [];
        columns_sarco_link = 53:56:56*this.nfibres;
    end
    
    % JLamDot
    if ~isempty(this.lambda_dot_pos)
        ildot = [];
        jldot = [];
    end

    %% -I part in u'(t) = -v(t)
    i = (1:3*N)';
    j = ((1:3*N)+3*N)';

    globidx_disp = sys.idx_u_glob_elems;
    globidx_press = sys.idx_p_glob_elems;

    dofs_pos = N*3;
    dofsperelem_displ = geo.DofsPerElement;
    dofsperelem_press = pgeo.DofsPerElement;
    num_gausspoints = fe_pos.GaussPointsPerElem;
    num_elements = geo.NumElements;
    pones = ones(dofsperelem_press,1,'int32');
    for m = 1:num_elements
        elemidx_u = globidx_disp(:,:,m);
        elemidx_v = elemidx_u + dofs_pos;
        elemidx_p = globidx_press(:,m);
        inew = elemidx_v(:);
        one = ones(size(inew),'int32');
        for gp = 1:num_gausspoints
            for k = 1:dofsperelem_displ
                %% Grad_u K(u,v,w)
                % xdim
                i = [i; inew]; %#ok<*AGROW>
                j = [j; one*elemidx_u(1,k)];
                % ydim
                i = [i; inew]; 
                j = [j; one*elemidx_u(2,k)]; 
                % zdim
                i = [i; inew]; 
                j = [j; one*elemidx_u(3,k)]; 

                %% Grad_v K(u,v,w)
%                 if visc > 0
%                     % xdim
%                     i = [i; inew]; %#ok<*AGROW>
%                     j = [j; one*elemidx_velo(1,k)];
%                     % ydim
%                     i = [i; inew]; 
%                     j = [j; one*elemidx_velo(2,k)]; 
%                     % zdim
%                     i = [i; inew]; 
%                     j = [j; one*elemidx_velo(3,k)]; 
%                 end

                %% grad u g(u)
                % dx
                i = [i; elemidx_p(:)];
                j = [j; pones*elemidx_u(1,k)]; 
                % dy
                i = [i; elemidx_p(:)];
                j = [j; pones*elemidx_u(2,k)]; 
                %dz
                i = [i; elemidx_p(:)];
                j = [j; pones*elemidx_u(3,k)]; 
            end
            %% Grad_w K(u,v,w)
            inew = elemidx_v(:);
            for k = 1:dofsperelem_press
                i = [i; inew];
                j = [j; ones(3*dofsperelem_displ,1,'int32')*elemidx_p(k)]; 
            end
            
            %% Jalpha pattern
            for k = 1:this.nfibres
                iS = [iS; elemidx_v(:)-dofs_pos];
                jS = [jS; ones(3*dofsperelem_displ,1)*columns_sarco_link(k)];
            end
            
            %% Check if change rate of lambda at a certain point should be tracked
            if ~isempty(this.lambda_dot_pos)
                k = find(this.lambda_dot_pos(1,:) == m & this.lambda_dot_pos(2,:) == gp);
                if ~isempty(k)
                    ildot = [ildot; k*ones(6*dofsperelem_displ,1)];
                    jldot = [jldot; elemidx_u(:); elemidx_v(:)];
                end
            end
        end
    end
    SP = sparse(double(i),double(j),ones(size(i)),6*N+M,6*N+M);
    % Remove values at dirichlet nodes
    SP(:,sys.idx_uv_bc_glob) = [];
    SP(sys.idx_uv_bc_glob,:) = [];
    SP = logical(SP);

    SPalpha = [];
    if this.nfibres > 0
        SPalpha = sparse(double(iS),double(jS),ones(size(iS)),3*N,this.nfibres*56);
        % Remove those that are connected to dirichlet values
        SPalpha([sys.idx_u_bc_glob; sys.idx_v_bc_glob],:) = [];
        SPalpha = logical(SPalpha);
    end
    
    SPLamDot = [];
    if ~isempty(this.lambda_dot_pos)
        SPLamDot = sparse(ildot,double(jldot),true(size(ildot)),this.nfibres,6*N);
        SPLamDot(:,sys.idx_uv_bc_glob) = [];
    end
end

