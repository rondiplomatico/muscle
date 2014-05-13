classdef System < models.BaseDynSystem
% MuscleFibreSystem: The global dynamical system used within the MuscleFibreModel
%
% Contains the FibreDynamics, Linear Diffusion and InputConv components as well as a
% ConstInitialValue.
%
% @author Daniel Wirtz @date 2014-01-20

    properties
       globidx_displ;
       
       globidx_pressure;
       
       Plota0 = true;
    end
    
    properties(SetAccess=private)
       pressure_to_displ_nodes;
       
       % The overall values of dirichlet conditions
       %
       % To specify in [mm] for position and mm/ms for velocity
       % 
       % Collected from bc_dir_displ_val, bc_dir_velo_val
       bc_dir_val;
       
       % The positions of the dirichlet values within the GLOBAL node/xyz
       % state space vector
       %
       % Collected from bc_dir_displ_idx, bc_dir_velo_idx
       bc_dir_idx;
       
       % Boundary conditions: 3 times numnodes logical matrix telling
       % which degree of freedom of which (position) node is fixed
       %     n1   n2   n3 ...
       % x | 1    0    0       |
       % y | 1    1    0 ...   |
       % z | 1    0    0       |
       bc_dir_displ;
       bc_dir_displ_idx;
       bc_dir_displ_val; % [mm]
       
       bc_dir_velo;
       bc_dir_velo_idx;
       bc_dir_velo_val; % [mm/ms]
       
       % A helper array containing the indices of the actual dofs in the
       % global indexing.
       %
       % Used to join dofs with dirichlet values
       dof_idx_global;
       
       % The indices of velocity components in the effective dof vector
       dof_idx_displ;
       dof_idx_velo;
       
       % Flag to invert the velocity mass matrix before simulations.
       %
       % If this is set to true, the mass matrix for the velocity part will
       % be assembled, the boundary condition rows removed and the inverse
       % of the remaining matrix will be pre-computed.
       % This inverse will be pre-multiplied to the velocity-dofs inside
       % the muscle.Dynamics evaluate and getStateJacobian functions.
       %
       % This in general results in higher simulation speed but it remains
       % open to see how well reduced modeling will work with that scheme.
       %
       % @type logical @default false
       UseDirectMassInversion = false;
       Minv;
       
       % Fibre stuff
       HasFibres = false;
       a0;
       dtna0;
       a0oa0;
       dNa0;
    end
    
    methods
        function this = System(model)
            % The system x'' + K(x) is transformed into the first order system
            % v' + K(x), x' = v
            
            % Call superclass constructor
            this = this@models.BaseDynSystem(model);
            
            this.addParam('mean input current',[0 1],10);
            
            %% Set system components
            % Core nonlinearity
            this.f = muscle.Dynamics(this);
        end
        
        function configUpdated(this)
            mc = this.Model.Config;
            tq = mc.PosFE;
            g = tq.Geometry;
            tl = mc.PressFE;
            gp = tl.Geometry;
                
%             % Find the indices of the pressure nodes in the displacement
%             % nodes geometry (used for plotting)
            this.pressure_to_displ_nodes = g.getCommonNodesWith(gp);
            
            % Call subroutine for boundary condition index crunching
            this.computeBC;
            
            % Init fibre directions and precomputable values
            this.inita0;
            
            % Construct global indices in uvw from element nodes. Each dof in
            % an element is used three times for x,y,z displacement. The
            % "elems" matrix contains the overall DOF numbers of each
            % element in the order of the nodes (along row) in the master
            % element.
            ne = g.NumElements;
            globalelementdofs = zeros(3,g.DofsPerElement,ne);
            for m = 1:ne
                % First index of element dof in global array
                hlp = (g.Elements(m,:)-1)*3+1;
                % Use first, second and third as positions.
                globalelementdofs(:,:,m) = [hlp; hlp+1; hlp+2];
            end
            this.globidx_displ = globalelementdofs;
            
            % The same for the pressure
            globalpressuredofs = zeros(gp.DofsPerElement,gp.NumElements);
            off = g.NumNodes * 6;
            for m = 1:gp.NumElements
                globalpressuredofs(:,m) = off + gp.Elements(m,:);
            end
            this.globidx_pressure = globalpressuredofs;
            
            %% Compile Mass Matrix
            % Augment mass matrix for all 3 displacement directions
            nd = g.NumNodes;
            [i, j, s] = find(tq.M);
            I = [3*(i'-1)+1; 3*(i'-1)+2; 3*(i'-1)+3];
            J = [3*(j'-1)+1; 3*(j'-1)+2; 3*(j'-1)+3];
            S = repmat(1:length(s),3,1);
            MM = sparse(I(:),J(:),s(S(:)),3*nd,3*nd);
            % Insert identity for velocity and zeros for pressure
            MM = blkdiag(speye(size(MM)),...
                MM,sparse(gp.NumNodes,gp.NumNodes));
            
            % Strip out the entries of dirichlet nodes
            MM(this.bc_dir_idx,:) = [];
            MM(:,this.bc_dir_idx) = [];
            
            % See description of property
            if this.UseDirectMassInversion
                this.Minv = inv(MM(this.dof_idx_velo,this.dof_idx_velo));
                MM = sparse(I(:),J(:),s(S(:)),3*nd,3*nd);
                % Use identity on left hand side
                MM = blkdiag(speye(size(MM)),...
                    speye(size(MM)),sparse(gp.NumNodes,gp.NumNodes));
                MM(this.bc_dir_idx,:) = [];
                MM(:,this.bc_dir_idx) = [];
            end
            this.M = dscomponents.ConstMassMatrix(MM);
            
            %% Initial value
            this.x0 = this.assembleX0;
            
            this.f.configUpdated;
        end
        
        function pm = plot(this, t, uvw, withvelo, pm, vid)
            if nargin < 6
                vid = false;
                if nargin < 4
                    withvelo = true;
                end
            end
            if nargin < 5 || isempty(pm)
                pm = PlotManager;
                pm.LeaveOpen = true;

            end            
            mc = this.Model.Config;
            dfem = mc.PosFE;
            geo = dfem.Geometry;
            pfem = mc.PressFE;
            pgeo = pfem.Geometry;
            vstart = geo.NumNodes * 3+1;
            pstart = geo.NumNodes * 6+1;
            e = geo.Edges;
            
            %% Re-add the dirichlet nodes
            yall = zeros(geo.NumNodes * 6 + pgeo.NumNodes, size(uvw,2));
            yall(this.dof_idx_global,:) = uvw;
            for k=1:length(this.bc_dir_idx)
                yall(this.bc_dir_idx(k),:) = this.bc_dir_val(k);
            end
            uvw = yall;
            
            bc_dir_displ_applies = sum(this.bc_dir_displ,1) >= 1;
            bc_dir_velo_applies = sum(this.bc_dir_velo,1) >= 1;
            
            if vid
                avifile = fullfile(pwd,'output.avi');
                vw = VideoWriter(avifile);
                vw.FrameRate = 30;
                vw.open;
            end
            
            %% Loop over time
            h = pm.nextPlot('geo','Output','x [mm]','y [mm]');
            view(h, [46 30]);
            daspect([1 1 1]);
            zlabel(h,'z [mm]');
            hold(h,'on');
%             box2 = dfem.geo.getBoundingBox(1.1);
            
            xpos = 1:3:geo.NumNodes*3;
            box = [min(min(uvw(xpos,:))) max(max(uvw(xpos,:)))...
                min(min(uvw(xpos+1,:))) max(max(uvw(xpos+1,:)))...
                min(min(uvw(xpos+2,:))) max(max(uvw(xpos+2,:)))]*1.1;
            for ts = 1:1:length(t)
                % Quit if figure has been closed
                if ~ishandle(h)
                    break;
                end
                u = reshape(uvw(1:vstart-1,ts),3,[]);
                v = reshape(uvw(vstart:pstart-1,ts),3,[]);
                cla(h);
                plot3(h,u(1,:),u(2,:),u(3,:),'k.','MarkerSize',14);
                
                % Velocities
                if withvelo
                    %quiver3(h,u(1,:),u(2,:),u(3,:),v(1,:),v(2,:),v(3,:),'k.','MarkerSize',14);
%                     quiver3(h,u(1,:),u(2,:),u(3,:),v(1,:),v(2,:),v(3,:),0,'r','MarkerSize',10);
                    quiver3(h,u(1,:),u(2,:),u(3,:),v(1,:),v(2,:),v(3,:),'b.', 'MarkerSize',14);
                end
                
                %% Dirichlet conditions
                % Displacement
                plot3(h,u(1,bc_dir_displ_applies),u(2,bc_dir_displ_applies),u(3,bc_dir_displ_applies),'k.','MarkerSize',20);
                for k=1:size(e,1)
                    plot3(h,u(1,[e(k,1) e(k,2)]),u(2,[e(k,1) e(k,2)]),u(3,[e(k,1) e(k,2)]),'r');
                end
                % Velocity
                plot3(h,u(1,bc_dir_velo_applies),u(2,bc_dir_velo_applies),u(3,bc_dir_velo_applies),'g.','MarkerSize',20);
                
                %% Pressure
                p = uvw(pstart:end,ts);
                pn = this.pressure_to_displ_nodes;
                pneg = p<0;
                % Negative pressures
                if any(pneg)
                    scatter3(h,u(1,pn(pneg)),u(2,pn(pneg)),u(3,pn(pneg)),-p(pneg),'r');
                end
                % Positive pressures
                if any(~pneg)
                    scatter3(h,u(1,pn(~pneg)),u(2,pn(~pneg)),u(3,pn(~pneg)),p(~pneg),'b');
                end
                
                %% a0 fibres
                if this.HasFibres && this.Plota0
                    Ngp = dfem.N(geo.gaussp);
                    for m = 1:geo.NumElements
                        u = uvw(1:vstart-1,ts);
                        u = u(this.globidx_displ(:,:,m));
                        gps = u*Ngp;
                        anull = u*this.dNa0(:,:,m);
                        quiver3(gps(1,:),gps(2,:),gps(3,:),anull(1,:),anull(2,:),anull(3,:),.5,'g');
                    end
                end
                
                %% Misc
                %axis(h,[-1.3 1.3 -1.3 1.3 -1.3 1.3]);
                axis(h,box);
%                view(h, [46 30]);
                title(h,sprintf('Deformation at t=%g',t(ts)));
%                 hold(h,'off');
                
                if vid
                    vw.writeVideo(getframe(gcf));
                else
%                     pause(.05);
                    pause;
                end
            end
            
            if vid
                vw.close;
            end

            if nargin < 4
                pm.done;
            end
        end
        
        function pm = plotDiff(this, t, uvw1, uvw2, fac, varargin)
            if nargin < 5
                fac = 5;
            end
            x0 = this.x0.evaluate([]);
            diff = repmat(x0,1,length(t)) + (uvw1-uvw2)*fac;
            pm = this.plot(t,diff,varargin{:});
        end
    end
    
    methods(Access=private)
        function x0 = assembleX0(this)
            % Constant initial values as current node positions
            mc = this.Model.Config;
            tq = mc.PosFE;
            geo = tq.Geometry;
            tl = mc.PressFE;
            pgeo = tl.Geometry;
            % All zero, especially the first tq.NumNodes * 3 velocity
            % entries and pressure
            x0 = zeros(geo.NumNodes * 6 + pgeo.NumNodes,1);
            
            % Fill in the reference configuration positions as initial
            % conditions
            for m = 1:geo.NumElements
                 dofpos = this.globidx_displ(:,:,m);
                 x0(dofpos) = geo.Nodes(:,geo.Elements(m,:));
            end
            
            % Initial conditions for pressure (s.t. S(X,0) = 0)
            x0(geo.NumNodes * 6+1:end) = -2*this.f.c10-4*this.f.c01;
            
%             velo_dir = false(3,tq.NumNodes);  
%             for k = [1 2 7 8]
%             for k = 1
%                 % Quadratic
%                 velo_dir(1,tq.elems(k,[1:3 9 10 13:15])) = true;
%                 %velo_dir(3,tq.elems(1,[1:3 9 10 13:15])) = true;
%     %             velo_dir(1,tq.elems(1,[3 10 15])) = true;
%             end
%             x0(find(velo_dir)+tq.NumNodes * 3) = .5;
            
            % Remove dirichlet values
            x0(this.bc_dir_idx) = [];
            
            x0 = dscomponents.ConstInitialValue(x0);
        end
        
        function computeBC(this)
            mc = this.Model.Config;
            [displ_dir, velo_dir, velo_dir_val] = mc.getBC;
            
            fe_displ = mc.PosFE;
            geo = fe_displ.Geometry;
            fe_press = mc.PressFE;
            pgeo = fe_press.Geometry;
            
            % Position of position entries in global state space vector
            num_displacement_dofs = geo.NumNodes * 3;
            
            %% Displacement
            this.bc_dir_displ = displ_dir;
            % Set values to node positions
            this.bc_dir_displ_val = geo.Nodes(displ_dir);
            % Add zeros for respective velocities
            this.bc_dir_displ_val = [this.bc_dir_displ_val; zeros(size(this.bc_dir_displ_val))];
            % Collect indices of dirichlet values per 3-set of x,y,z values
            relpos = find(displ_dir(:));
            % Same positions for points and velocity
            this.bc_dir_displ_idx = [relpos; num_displacement_dofs + relpos];
            
            %% Velocity
            this.bc_dir_velo = velo_dir;
            this.bc_dir_velo_val = velo_dir_val(velo_dir);
            this.bc_dir_velo_idx = num_displacement_dofs + find(velo_dir(:));
            
            %% Hydrostatic Pressure Dirichlet conditions
%             press_dir = false(1,tl.NumNodes);
%             press_dir(1) = true;
%             this.bc_dir = [displ_dir; press_dir];

            % The according velocities are all zero for dirichlet points
            this.bc_dir_val = [this.bc_dir_displ_val; this.bc_dir_velo_val];
            this.bc_dir_idx = [this.bc_dir_displ_idx; this.bc_dir_velo_idx];
            
            % Compute dof positions in global state space vector
            total = geo.NumNodes * 6 + pgeo.NumNodes;
            pos = false(1,total);
            pos(1:num_displacement_dofs) = true;
            pos(this.bc_dir_idx) = [];
            this.dof_idx_displ = find(pos);
            
            pos = false(1,total);
            pos(num_displacement_dofs+1:num_displacement_dofs*2) = true;
            pos(this.bc_dir_idx) = [];
            this.dof_idx_velo = find(pos);
            
            idx = 1:total;
            idx(this.bc_dir_idx) = [];
            this.dof_idx_global = idx;
        end
        
        function inita0(this)
            mc = this.Model.Config;
            fe = mc.PosFE;
            geo = fe.Geometry;
            
            anull = mc.geta0;
            this.HasFibres = false;
            if any(anull(:))
                this.HasFibres = true;
            
                this.a0 = anull;

                % Precomputations
                dNgp = fe.gradN(geo.gaussp);
                anulldyadanull = zeros(3,3,geo.GaussPointsPerElem*geo.NumElements);
                dtnanull = zeros(geo.DofsPerElement,geo.GaussPointsPerElem,geo.NumElements);
                dNanull = zeros(geo.DofsPerElement,geo.GaussPointsPerElem,geo.NumElements);
                for m = 1 : geo.NumElements
                    for gp = 1 : geo.GaussPointsPerElem
                        % a0 dyad a0
                        pos = (m-1)*geo.NumElements+gp;
                        anulldyadanull(:,:,pos) = anull(:,gp,m)*anull(:,gp,m)';

                        % <grad phi_k, a0> scalar products
                        pos = 3*(gp-1)+1:3*gp;
                        dtn = fe.transgrad(:,pos,m);
                        dtnanull(:,gp,m) = dtn*anull(:,gp,m);

                        % forward transformation of a0 at gauss points
                        % (plotting only so far)
                        pos = [0 27 54]+gp;
                        dNanull(:,gp,m) = dNgp(:,pos) * this.a0(:,gp,m);
                    end
                end
                this.dtna0 = dtnanull;
                this.a0oa0 = anulldyadanull;
                this.dNa0 = dNanull;
            end
        end
    end
    
end
