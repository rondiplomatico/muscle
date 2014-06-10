classdef System < models.BaseDynSystem
% MuscleFibreSystem: The global dynamical system used within the MuscleFibreModel
%
% The system u'' + K(u) is transformed into the first order system v' +
% K(u), u' = v
%
% @author Daniel Wirtz @date 2014-01-20

    properties
       globidx_displ;
       
       globidx_pressure;
       
       % Set this to a double value to apply velocity dirichlet conditions
       % only up to a certain time (zero after that)
       ApplyVelocityBCUntil = Inf;
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
              
       bc_neum_forces_nodeidx; % [N]
       bc_neum_forces_val;
       FacesWithForce;
       
       % A helper array containing the indices of the actual dofs in the
       % global indexing.
       %
       % Used to join dofs with dirichlet values
       dof_idx_global;
       
       % The indices of velocity components in the effective dof vector
       dof_idx_displ;
       dof_idx_velo;
       
       Minv;
       
       % Fibre stuff
       HasFibres = false;
       a0;
       dtna0;
       a0oa0;
       dNa0;
    end
   
    properties(Dependent)
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
       UseDirectMassInversion;
    end
    
    properties(Access=private)
        fUseDirectMassInversion = false;
        fD;
    end
    
    methods
        function this = System(model)
            % Call superclass constructor
            this = this@models.BaseDynSystem(model);
            
            this.addParam('viscosity',[0 10],10);
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
            this.computeDirichletBC;
            
            %% Construct B matrix
            % Collect neumann forces
            [B, this.bc_neum_forces_nodeidx] = this.getSpatialExternalForces;
            % Only set up forces if present
            if ~isempty(this.bc_neum_forces_nodeidx)
                this.bc_neum_forces_val = B(this.bc_neum_forces_nodeidx + g.NumNodes * 3);
                % Remove dirichlet DoFs
                B(this.bc_dir_idx) = [];
                % Set as constant input conversion matrix
                this.B = dscomponents.LinearInputConv(B);
                % Set input function
                this.Inputs{1} = mc.getInputFunction(this.Model);
            end
            
            % Init fibre directions and precomputable values
            this.inita0;
            
            % Construct global indices in uvw from element nodes. Each dof in
            % an element is used three times for x,y,z displacement. The
            % "elems" matrix contains the overall DOF numbers of each
            % element in the order of the nodes (along row) in the master
            % element.
            ne = g.NumElements;
            globalelementdofs = zeros(3,g.DofsPerElement,ne,'int32');
            for m = 1:ne
                % First index of element dof in global array
                hlp = (g.Elements(m,:)-1)*3+1;
                % Use first, second and third as positions.
                globalelementdofs(:,:,m) = [hlp; hlp+1; hlp+2];
            end
            this.globidx_displ = globalelementdofs;
            
            % The same for the pressure
            globalpressuredofs = zeros(gp.DofsPerElement,gp.NumElements,'int32');
            off = g.NumNodes * 6;
            for m = 1:gp.NumElements
                globalpressuredofs(:,m) = off + gp.Elements(m,:);
            end
            this.globidx_pressure = globalpressuredofs;
            
            %% Compile Mass Matrix
            this.M = this.assembleMassMatrix;
            
            %% Compile Damping Matrix
            this.fD = this.assembleDampingMatrix;
            
            %% Initial value
            this.x0 = this.assembleX0;
            
            this.f.configUpdated;
        end
        
        function prepareSimulation(this, mu, inputidx)
            this.A = [];
            if mu(1) > 0
                this.A = this.fD;
            end
            prepareSimulation@models.BaseDynSystem(this, mu, inputidx);
        end
        
        function pm = plot(this, t, uvw, varargin)
            i = inputParser;
            i.KeepUnmatched = true;
            i.addParamValue('Vid',false,@(v)islogical(v));
            i.addParamValue('Forces',false,@(v)islogical(v));
            i.addParamValue('Velo',false,@(v)islogical(v));
            i.addParamValue('Pressure',false,@(v)islogical(v));
            i.addParamValue('Fibres',true,@(v)islogical(v));
            i.addParamValue('Skel',false,@(v)islogical(v));
            i.addParamValue('Pool',length(t)>1,@(v)islogical(v));
            i.addParamValue('pm',[],@(v)isa(v,'PlotManager'));
            i.addParamValue('DF',[]);
            i.addParamValue('NF',[]);
            i.parse(varargin{:});
            r = i.Results;
            if ~isempty(r.NF)
                r.Forces = true;
            end
            
            mc = this.Model.Config;
            
            if isempty(r.pm)
                if ~isempty(mc.Pool) && r.Pool
                    pm = PlotManager(false,2,1);
                else
                    pm = PlotManager;
                end
                pm.LeaveOpen = true;
            else
                pm = r.pm;
            end
            
            dfem = mc.PosFE;
            geo = dfem.Geometry;
            posdofs = geo.NumNodes * 3;
            vstart = posdofs+1;
            pstart = 2*posdofs+1;
            e = geo.Edges;
            
            %% Re-add the dirichlet nodes
            uvw = this.includeDirichletValues(t, uvw);
            
            hlp = sum(this.bc_dir_displ,1);
            bc_dir_3pos_applies = hlp == 3;
            bc_dir_2pos_applies = hlp == 2;
            bc_dir_1pos_applies = hlp == 1;
            bc_dir_pos_applies = hlp >= 1; 
            bc_dir_velo_applies = sum(this.bc_dir_velo,1) >= 1;
            no_bc = ~bc_dir_pos_applies & ~bc_dir_velo_applies;
            
            if ~isempty(r.DF)
                have_residuals = bc_dir_pos_applies | bc_dir_velo_applies;
                % By sorting and combining the pos/velo Dir BC, the
                % plotting of mixed BCs on one node is plotted correctly.
                residuals_pos = this.bc_dir_displ | this.bc_dir_velo;
                [~, sortidx] = sort([this.bc_dir_displ_idx; this.bc_dir_velo_idx-posdofs]);
                % Preallocate the residuals matrix
                residuals = zeros(size(residuals_pos));
                maxdfval = max(abs(r.DF(:)))/10;
            end
            
            if r.Vid
                avifile = fullfile(pwd,'output.avi');
                vw = VideoWriter(avifile);
                vw.FrameRate = 30;
                vw.open;
            end
            
            if r.Forces
                % Get forces on each node in x,y,z directions
                % This is where the connection between plane index and
                % x,y,z coordinate is "restored"
                forces = zeros(size(this.bc_dir_displ));
                forces(this.bc_neum_forces_nodeidx) = this.bc_neum_forces_val;
                
                % Forces at face centers
                force_elem_face_idx = geo.Faces(:,this.FacesWithForce);
                numfaceswithforce = size(force_elem_face_idx,2);
                meanforces = zeros(3,numfaceswithforce);
                for k=1:numfaceswithforce
                    masterfacenodeidx = geo.MasterFaces(force_elem_face_idx(2,k),:);
                    facenodeidx = geo.Elements(force_elem_face_idx(1,k),masterfacenodeidx);
                    meanforces(:,k) = mean(forces(:,facenodeidx),2);
                end
                
                % Also save forces at x,y,z nodes for plotting
                forces_apply = sum(abs(forces),1) ~= 0;
                forces = forces(:,forces_apply);
                
                if ~isempty(r.NF)
                    residual_neumann_forces = zeros(size(this.bc_dir_displ));
                end
            end
            
            %% Loop over time
            if ~isempty(mc.Pool) && r.Pool
                pool = mc.Pool;
                ha = pm.nextPlot('force','Activation force','t [ms]','alpha');
                axis(ha,[0 t(end) 0 1]);
                hold(ha,'on');
                dt = this.Model.dt;
                typeweights = mc.FibreTypeWeights(1,:,1);
            end
            
            h = pm.nextPlot('geo','Output','x [mm]','y [mm]');
            zlabel(h,'z [mm]');
            
            if ~r.Skel
%                 light('Position',[1 1 1],'Style','infinite','Parent',h);
                musclecol = [0.854688, 0.201563, 0.217188];
            end
            
            axis(h, this.getPlotBox(uvw));
            daspect([1 1 1]);
            view(h, [46 30]);
            hold(h,'on');
            
            for ts = 1:length(t)
                % Quit if figure has been closed
                if ~ishandle(h)
                    break;
                end
                u = reshape(uvw(1:vstart-1,ts),3,[]);
                v = reshape(uvw(vstart:pstart-1,ts),3,[]);
                cla(h);
                
                if r.Skel
                    plot3(h,u(1,no_bc),u(2,no_bc),u(3,no_bc),'r.','MarkerSize',14);
                    for k=1:size(e,1)
                        plot3(h,u(1,[e(k,1) e(k,2)]),u(2,[e(k,1) e(k,2)]),u(3,[e(k,1) e(k,2)]),'r');
                    end
                else
                    p = patch('Faces',geo.PatchFaces,'Vertices',u');
                    set(p,'EdgeColor',.8*musclecol,'FaceColor',musclecol,'FaceAlpha',.3);
                end
                
                % Velocities
                if r.Velo
                    quiver3(h,u(1,:),u(2,:),u(3,:),v(1,:),v(2,:),v(3,:),'b', 'MarkerSize',14);
                end
                
                %% Dirichlet conditions
                % Displacement
                plot3(h,u(1,bc_dir_3pos_applies),u(2,bc_dir_3pos_applies),u(3,bc_dir_3pos_applies),'.','MarkerSize',20,'Color',[0 0 0]);
                plot3(h,u(1,bc_dir_2pos_applies),u(2,bc_dir_2pos_applies),u(3,bc_dir_2pos_applies),'.','MarkerSize',20,'Color',[.5 .5 .5]);
                plot3(h,u(1,bc_dir_1pos_applies),u(2,bc_dir_1pos_applies),u(3,bc_dir_1pos_applies),'.','MarkerSize',20,'Color',[.7 .7 .7]);
                
                % Velocity
                plot3(h,u(1,bc_dir_velo_applies),u(2,bc_dir_velo_applies),u(3,bc_dir_velo_applies),'g.','MarkerSize',20);
                
                %% Dirichlet Forces
                if ~isempty(r.DF)
                    udir = u(:,have_residuals);
                    residuals(residuals_pos) = r.DF(sortidx,ts)/maxdfval;
                    quiver3(h,udir(1,:),udir(2,:),udir(3,:),...
                        residuals(1,have_residuals),residuals(2,have_residuals),residuals(3,have_residuals),'k', 'MarkerSize',14);
                end
                
                %% Neumann condition forces
                if r.Forces
                    % Plot force vectors at nodes
                    uforce = u(:,forces_apply);
                    quiver3(h,uforce(1,:),uforce(2,:),uforce(3,:),...
                        forces(1,:),forces(2,:),forces(3,:),0,'Color',[.8 .8 1]);
                    % Plot force vectors at face centers
                    for k=1:numfaceswithforce
                        masterfacenodeidx = geo.MasterFaces(force_elem_face_idx(2,k),:);
                        facenodeidx = geo.Elements(force_elem_face_idx(1,k),masterfacenodeidx);
                        
                        facecenter = mean(u(:,facenodeidx),2);
                        quiver3(h,facecenter(1),facecenter(2),facecenter(3),...
                            meanforces(1,k),meanforces(2,k),meanforces(3,k),'LineWidth',2,'Color','b','MaxHeadSize',1);
                        
                        if ~isempty(r.NF)
                            residual_neumann_forces(this.bc_neum_forces_nodeidx) = r.NF(:,ts);
                            meanforce = mean(residual_neumann_forces(:,facenodeidx),2);
                            quiver3(h,facecenter(1),facecenter(2),facecenter(3),...
                            meanforce(1),meanforce(2),meanforce(3),'LineWidth',2,'Color','k','MaxHeadSize',1);
                        end
                    end
                end
                
                %% Pressure
                if r.Pressure
                    p = uvw(pstart:end,ts);
                    pn = this.pressure_to_displ_nodes;
                    pneg = p<0;
                    % Negative pressures
                    if any(pneg)
                        scatter3(h,u(1,pn(pneg)),u(2,pn(pneg)),u(3,pn(pneg)),-p(pneg)*10,'b');
                    end
                    % Positive pressures
                    if any(~pneg)
                        scatter3(h,u(1,pn(~pneg)),u(2,pn(~pneg)),u(3,pn(~pneg)),p(~pneg)*10,'b');
                    end
                end
                
                %% a0 fibres
                if this.HasFibres && r.Fibres
                    Ngp = dfem.N(dfem.GaussPoints);
                    for m = 1:geo.NumElements
                        u = uvw(1:vstart-1,ts);
                        u = u(this.globidx_displ(:,:,m));
                        
                        gps = u*Ngp;
                        anull = u*this.dNa0(:,:,m);
                        quiver3(gps(1,:),gps(2,:),gps(3,:),anull(1,:),anull(2,:),anull(3,:),.5,'.','Color','w');
                    end
                end
                
                %% Misc
%                 axis(h,box);
%                view(h, [73 40]);
%                 view(h, [0 90]);
%                 view(h, [-90 0]);
                title(h,sprintf('Deformation at t=%g',t(ts)));
%                 hold(h,'off');

                if ~isempty(mc.Pool) && r.Pool
                    times = 0:dt:ts*dt;
                    alpha = pool.getActivation(times);
                    walpha = typeweights * alpha;
                    cla(ha);
                    plot(ha,times,alpha);
                    plot(ha,times,walpha,'LineWidth',2);
                end
                
                if r.Vid
                    vw.writeVideo(getframe(gcf));
                else
                    pause(.05);
%                     pause;
                end
            end
            
            if r.Vid
                vw.close;
            end

            if isempty(r.pm)
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
        
        function uvwall = includeDirichletValues(this, t, uvw)
            mc = this.Model.Config;
            dfem = mc.PosFE;
            geo = dfem.Geometry;
            pfem = mc.PressFE;
            pgeo = pfem.Geometry;
            
            %% Re-add the dirichlet nodes
            uvwall = zeros(geo.NumNodes * 6 + pgeo.NumNodes, size(uvw,2));
            uvwall(this.dof_idx_global,:) = uvw;
            uvwall(this.bc_dir_idx,:) = repmat(this.bc_dir_val,1,size(uvw,2));
            
            sys = this.Model.System;
            zerovel = t > sys.ApplyVelocityBCUntil;
            uvwall(sys.bc_dir_velo_idx,zerovel) = 0;
        end
        
    end
    
    methods
        function set.UseDirectMassInversion(this, value)
            if ~islogical(value) || ~isscalar(value)
                error('UseDirectMassInversion must be true or false');
            end
            if this.fUseDirectMassInversion ~= value
                this.fUseDirectMassInversion = value;
                this.configUpdated;
            end
        end
        
        function value = get.UseDirectMassInversion(this)
            value = this.fUseDirectMassInversion;
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

            % Give the model config a chance to provide x0
            x0 = mc.getX0(x0);

            % Remove dirichlet values
            x0(this.bc_dir_idx) = [];
            
            x0 = dscomponents.ConstInitialValue(x0);
        end
        
        function M = assembleMassMatrix(this)
            %% Compile Mass Matrix
            
            mc = this.Model.Config;
            fe_pos = mc.PosFE;
            g = fe_pos.Geometry;
            gp = mc.PressFE.Geometry;
            
            % Augment mass matrix for all 3 displacement directions
            nd = g.NumNodes;
            [i, j, s] = find(fe_pos.M);
            I = [3*(i'-1)+1; 3*(i'-1)+2; 3*(i'-1)+3];
            J = [3*(j'-1)+1; 3*(j'-1)+2; 3*(j'-1)+3];
            S = repmat(1:length(s),3,1);
            MM = sparse(I(:),J(:),s(S(:)),3*nd,3*nd);
            % Insert identity for velocity and zeros for pressure
            MM = blkdiag(speye(size(MM)),...
                this.Model.MuscleDensity*MM,sparse(gp.NumNodes,gp.NumNodes));
            
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
            M = dscomponents.ConstMassMatrix(MM);
        end
        
        function Daff = assembleDampingMatrix(this)
            %% Compile Damping/Viscosity Matrix
            mc = this.Model.Config;
            fe_pos = mc.PosFE;
            g = fe_pos.Geometry;
            gp = mc.PressFE.Geometry;
            
            % Augment mass matrix for all 3 displacement directions
            nd = g.NumNodes;
            [i, j, s] = find(fe_pos.D);
            I = [3*(i'-1)+1; 3*(i'-1)+2; 3*(i'-1)+3];
            J = [3*(j'-1)+1; 3*(j'-1)+2; 3*(j'-1)+3];
            S = repmat(1:length(s),3,1);
            D = sparse(I(:),J(:),s(S(:)),3*nd,3*nd);
            % Insert identity for position and zeros for pressure
            D = blkdiag(sparse(3*nd,3*nd),D,sparse(gp.NumNodes,gp.NumNodes));
            
            % Strip out the entries of dirichlet nodes
            D(this.bc_dir_idx,:) = [];
            D(:,this.bc_dir_idx) = [];
            
            Daff = dscomponents.AffLinCoreFun(this);
            Daff.addMatrix('mu(1)',-D);
        end
        
        function computeDirichletBC(this)
            mc = this.Model.Config;
            [pos_dir, velo_dir, velo_dir_val] = mc.getBC;
            
            fe_displ = mc.PosFE;
            geo = fe_displ.Geometry;
            fe_press = mc.PressFE;
            pgeo = fe_press.Geometry;
            
            % Position of position entries in global state space vector
            num_position_dofs = geo.NumNodes * 3;
            
            %% Displacement
            this.bc_dir_displ = pos_dir;
            % Set values to node positions
            this.bc_dir_displ_val = geo.Nodes(pos_dir);
            this.bc_dir_displ_idx = int32(find(pos_dir(:)));
            
            %% Velocity
            % Add any user-defines values (cannot conflict with position
            % dirichlet conditions, this is checked in AModelConfig.getBC)
            this.bc_dir_velo = velo_dir;
            this.bc_dir_velo_idx = int32(num_position_dofs + find(velo_dir(:)));
            this.bc_dir_velo_val = velo_dir_val(velo_dir);
            
            %% Hydrostatic Pressure Dirichlet conditions
%             press_dir = false(1,tl.NumNodes);
%             press_dir(1) = true;
%             this.bc_dir = [displ_dir; press_dir];

            % Compile the global dirichlet values index and value vectors.
            % Here we add zero velocities for each point with fixed position, too.
            [this.bc_dir_idx, sortidx] = sort([this.bc_dir_displ_idx; this.bc_dir_displ_idx+num_position_dofs; this.bc_dir_velo_idx]);
            alldirvals = [this.bc_dir_displ_val; zeros(size(this.bc_dir_displ_val)); this.bc_dir_velo_val];
            this.bc_dir_val = alldirvals(sortidx);
            
            % Compute dof positions in global state space vector
            total = geo.NumNodes * 6 + pgeo.NumNodes;
            pos = false(1,total);
            pos(1:num_position_dofs) = true;
            pos(this.bc_dir_idx) = [];
            this.dof_idx_displ = int32(find(pos));
            
            pos = false(1,total);
            pos(num_position_dofs+1:num_position_dofs*2) = true;
            pos(this.bc_dir_idx) = [];
            this.dof_idx_velo = int32(find(pos));
            
            idx = int32(1:total);
            idx(this.bc_dir_idx) = [];
            this.dof_idx_global = idx;
        end
        
        function [force, nodeidx] = getSpatialExternalForces(this)
            mc = this.Model.Config;
            fe_displ = mc.PosFE;
            geo = fe_displ.Geometry;
            ngp = fe_displ.GaussPointsPerElemFace;
            force = zeros(geo.NumNodes * 3,1);
            faceswithforce = false(1,geo.NumFaces);
            for fn = 1:geo.NumFaces
                elemidx = geo.Faces(1,fn);
                faceidx = geo.Faces(2,fn);
                masterfacenodeidx = geo.MasterFaces(faceidx,:);
                % So far: Constant pressure on all gauss points!
                P = mc.getBoundaryPressure(elemidx, faceidx);
                if ~isempty(P)
                    faceswithforce(fn) = true;
                    integrand = zeros(3,geo.NodesPerFace);
                    for gi = 1:ngp
                        PN = (P * fe_displ.NormalsOnFaceGP(:,gi,fn)) * fe_displ.Ngpface(:,gi,fn)';
                        integrand = integrand + fe_displ.FaceGaussWeights(gi)*PN*fe_displ.face_detjac(fn,gi);
                    end
                    facenodeidx = geo.Elements(elemidx,masterfacenodeidx);
                    facenodeidx = (facenodeidx-1)*3+1;
                    facenodeidx = [facenodeidx; facenodeidx+1; facenodeidx+2];%#ok
                    force(facenodeidx(:)) = force(facenodeidx(:)) + integrand(:);
                end
            end
            % Augment to u,v,w vector
            nodeidx = find(abs(force) > 1e-13);
            force = [zeros(size(force)); force; zeros(mc.PressFE.Geometry.NumNodes,1)];
            this.FacesWithForce = faceswithforce;
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
                dNgp = fe.gradN(fe.GaussPoints);
                
                anulldyadanull = zeros(3,3,fe.GaussPointsPerElem*geo.NumElements);
                dtnanull = zeros(geo.DofsPerElement,fe.GaussPointsPerElem,geo.NumElements);
                dNanull = zeros(geo.DofsPerElement,fe.GaussPointsPerElem,geo.NumElements);
                for m = 1 : geo.NumElements
                    u = geo.Nodes(:,geo.Elements(m,:));
                    for gp = 1 : fe.GaussPointsPerElem
                        
                        % Transform a0 to local fibre direction
                        pos = [0 fe.GaussPointsPerElem 2*fe.GaussPointsPerElem]+gp;
                        Jac = u*dNgp(:,pos);
                        loc_anull = Jac * anull(:,gp,m);
                        loc_anull = loc_anull/norm(loc_anull);
                        
                        % forward transformation of a0 at gauss points
                        % (plotting only so far)
                        dNanull(:,gp,m) = dNgp(:,pos) * anull(:,gp,m);
                        
                        % <grad phi_k, a0> scalar products (for
                        % getStateJacobian)
                        pos = 3*(gp-1)+1:3*gp;
                        dtnanull(:,gp,m) = fe.transgrad(:,pos,m)*loc_anull;
                        
                        % a0 dyad a0
                        pos = (m-1)*fe.GaussPointsPerElem+gp;
                        anulldyadanull(:,:,pos) = loc_anull*loc_anull';%anull(:,gp,m)*anull(:,gp,m)';
                    end
                end
                this.dtna0 = dtnanull;
                this.a0oa0 = anulldyadanull;
                this.dNa0 = dNanull;
            end
        end
        
        function box = getPlotBox(this, uvw)
            xpos = 1:3:this.Model.Config.PosFE.Geometry.NumNodes*3;
            box = [min(min(uvw(xpos,:))) max(max(uvw(xpos,:)))...
                min(min(uvw(xpos+1,:))) max(max(uvw(xpos+1,:)))...
                min(min(uvw(xpos+2,:))) max(max(uvw(xpos+2,:)))];
            diam = diff(box)*.1;
            box = box + [-diam(1) diam(1) -diam(3) diam(3) -diam(5) diam(5)];
        end
    end
    
end
