classdef System < models.BaseDynSystem
% MuscleFibreSystem: The global dynamical system used within the MuscleFibreModel
%
% The system u'' + K(u) is transformed into the first order system v' +
% K(u), u' = v
%
% @author Daniel Wirtz @date 2014-01-20

    properties
       % The global index of node x,y,z positions within uvw.
       %
       % The velocities are indentically indexed but 3*NumNodes later.
       idx_u_glob_elems;
       
       % The global index of node pressures within uvw.
       idx_p_glob_elems;
       
       % Set this to a double value to apply velocity dirichlet conditions
       % only up to a certain time (zero after that)
       %
       % This can be used to generate initial conditions with a deformed
       % state
       %
       % @type double @default Inf
       ApplyVelocityBCUntil = Inf;
    end
    
    properties(SetAccess=private)
       % This contains the indices of the nodes tracking pressure (linear /
       % 8-corner elements) within the quadratic global node numbering
       %
       % Used for plotting only so far.
       idx_p_to_u_nodes;
       
       % The overall values of dirichlet conditions
       %
       % To specify in [mm] for position and [mm/ms] for velocity
       % 
       % Collected from val_u_bc, val_v_bc
       val_uv_bc_glob;
       
       % The positions of the dirichlet values within the global node/xyz
       % state space vector
       %
       % Collected from idx_u_bc_glob, idx_v_bc_glob
       idx_uv_bc_glob;
       
       % Boundary conditions: 3 times numnodes logical matrix telling
       % which degree of freedom of which (position) node is fixed
       %     n1   n2   n3 ...
       % x | 1    0    0       |
       % y | 1    1    0 ...   |
       % z | 1    0    0       |
       bool_u_bc_nodes;
       idx_u_bc_glob;
       val_u_bc; % [mm]
       
       % Same
       bool_v_bc_nodes;
       idx_v_bc_glob;
       val_v_bc; % [mm/ms]
              
       bc_neum_forces_nodeidx; % [N]
       bc_neum_forces_val;
       FacesWithForce;
       
       % A helper array containing the indices of the actual dofs in the
       % global indexing.
       %
       % Used to join dofs with dirichlet values
       idx_uv_dof_glob;
       
       % The indices of velocity components in the effective dof vector
       idx_u_dof_glob;
       num_u_dof;
       idx_v_dof_glob;
       num_v_dof;
       idx_p_dof_glob;
       num_p_dof;
       
       % Total number of discrete field points including dirichlet values
       num_uvp_glob;
       % Total number of discrete field points that are degrees of freedom
       num_uvp_dof;
       
       Minv;
       
       %% Fibre stuff
       HasFibres = false;
       HasFibreTypes = false;
       HasMotoPool = false;
       % Property only useful when fullmodel.System is used.. comes due to
       % use of inheritance. Better solutions could be thought of, but this
       % is to get it going!
       HasForceArgument = false;
       a0;
       a0oa0;
       dNa0;
       a0Base;
       
       %% Cross-fibre stiffness stuff
       % normals to a0
       a0oa0n1;
       a0oa0n2;
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
       
       % Flag to indicate if this system should use stiffness terms
       % (Markert) for cross-fibre directions.
       %
       % @type logaical @default false
       UseCrossFibreStiffness;
    end
    
    properties(Access=private)
        fUseDirectMassInversion = false;
        fUseCrossFibreStiffness = false;
        fD;
    end
    
    methods
        function this = System(model)
            % Call superclass constructor
            this = this@models.BaseDynSystem(model);
            
            % The muscle viscosity. [mN * mm/ms]
            this.addParam('viscosity',[0 10],10);
            
            % The amount of milliseconds over which to activate the model
            this.addParam('alpha_ramp',[10,300],10);
            
            % The force of the faces exposed to neumann conditions
            this.addParam('neumann_force',[0 400],10);
            
            % For some variants, we have the mean input current for the
            % motoneuron pool (generating activation)
            this.addParam('mean input current',[0 1],10);
            
            %% Set system components
            % Core nonlinearity
            this.f = muscle.Dynamics(this);
        end
        
        function rsys = buildReducedSystem(this, rmodel)
            % Overrides the default buildReducedModel method in order to
            % temporarily set the A component used in projection.
            %
            % (So far the extra efficiency of "no A" for mu(1) == 0 is
            % neglected in any reduced model)
            this.A = this.fD;
            rsys = buildReducedSystem@models.BaseDynSystem(this, rmodel);
            this.A = [];
        end
        
        function configUpdated(this)
            mc = this.Model.Config;
            if ~isempty(mc)
                tq = mc.PosFE;
                geo_uv = tq.Geometry;
                tl = mc.PressFE;
                geo_p = tl.Geometry;

                % Find the indices of the pressure nodes in the displacement
                % nodes geometry (used for plotting)
                this.idx_p_to_u_nodes = geo_uv.getCommonNodesWith(geo_p);

                % Call subroutine for boundary condition index crunching
                this.computeDirichletBC;
                
                this.updateDofNums(mc);

                %% Construct B matrix
                % Collect neumann forces
                [B, this.bc_neum_forces_nodeidx] = this.getSpatialExternalForces;
                % Only set up forces if present
                if ~isempty(this.bc_neum_forces_nodeidx)
                    this.bc_neum_forces_val = B(this.bc_neum_forces_nodeidx + geo_uv.NumNodes * 3);
                    % Remove dirichlet DoFs
                    B(this.idx_uv_bc_glob) = [];
                    % Set as constant input conversion matrix
                    this.B = dscomponents.LinearInputConv(B);
                    % Set input function
                    this.Inputs = mc.getInputs;
                end

                % Init fibre directions and precomputable values
                this.inita0;
                
                this.HasMotoPool = this.HasFibres && ~isempty(mc.Pool);
                this.HasFibreTypes = this.HasFibres && ~isempty(mc.FibreTypeWeights);
                this.HasForceArgument = this.HasFibreTypes && isa(this,'fullmuscle.System');

                % Construct global indices in uvw from element nodes. Each dof in
                % an element is used three times for x,y,z displacement. The
                % "elems" matrix contains the overall DOF numbers of each
                % element in the order of the nodes (along row) in the master
                % element.
                ne = geo_uv.NumElements;
                globalelementdofs = zeros(3,geo_uv.DofsPerElement,ne,'int32');
                for m = 1:ne
                    % First index of element dof in global array
                    hlp = (geo_uv.Elements(m,:)-1)*3+1;
                    % Use first, second and third as positions.
                    globalelementdofs(:,:,m) = [hlp; hlp+1; hlp+2];
                end
                this.idx_u_glob_elems = globalelementdofs;

                % The same for the pressure
                globalpressuredofs = zeros(geo_p.DofsPerElement,geo_p.NumElements,'int32');
                off = geo_uv.NumNodes * 6;
                for m = 1:geo_p.NumElements
                    globalpressuredofs(:,m) = off + geo_p.Elements(m,:);
                end
                this.idx_p_glob_elems = globalpressuredofs;

                %% Compile Mass Matrix
                this.M = dscomponents.ConstMassMatrix(this.assembleMassMatrix);

                %% Compile Damping Matrix
                this.fD = this.assembleDampingMatrix;

                %% Initial value
                this.x0 = dscomponents.ConstInitialValue(this.assembleX0);

                this.f.configUpdated;
            end
        end
        
        function prepareSimulation(this, mu, inputidx)
            this.A = [];
            if mu(1) > 0
                this.A = this.fD;
            end
            prepareSimulation@models.BaseDynSystem(this, mu, inputidx);
        end
        
        function [pm, h] = plot(this, t, y_dofs, varargin)
            i = inputParser;
            i.KeepUnmatched = true;
            i.addParamValue('Vid',false,@(v)islogical(v));
            i.addParamValue('Forces',false,@(v)islogical(v));
            i.addParamValue('Velo',false,@(v)islogical(v));
            i.addParamValue('Pressure',false,@(v)islogical(v));
            i.addParamValue('Fibres',true,@(v)islogical(v));
            i.addParamValue('Skel',false,@(v)islogical(v));
            i.addParamValue('Pool',length(t)>1,@(v)islogical(v));
            i.addParamValue('PM',[],@(v)isa(v,'PlotManager'));
            i.addParamValue('DF',[]);
            i.addParamValue('NF',[]);
            i.addParamValue('F',[]);
            i.parse(varargin{:});
            r = i.Results;
            if ~isempty(r.NF)
                r.Forces = true;
            end
            
            mc = this.Model.Config;
            
            if isempty(r.PM)
                if ~isempty(mc.Pool) && r.Pool
                    pm = PlotManager(false,2,1);
                else
                    pm = PlotManager;
                end
                pm.LeaveOpen = true;
            else
                pm = r.PM;
            end
            
            dfem = mc.PosFE;
            geo = dfem.Geometry;
            posdofs = geo.NumNodes * 3;
            vstart = posdofs+1;
            pstart = 2*posdofs+1;
            e = geo.Edges;
            
            % "Speedup" factor for faster plots
            if ~isempty(r.F)
                t = t(1:r.F:end);
                y_dofs = y_dofs(:,1:r.F:end);
                if ~isempty(r.DF)
                    r.DF = r.DF(:,1:r.F:end);
                end
                if ~isempty(r.NF)
                    r.NF = r.NF(:,1:r.F:end);
                end
            end
            
            %% Re-add the dirichlet nodes
            y_dofs = this.includeDirichletValues(t, y_dofs);
            
            hlp = sum(this.bool_u_bc_nodes,1);
            bc_dir_3pos_applies = hlp == 3;
            bc_dir_2pos_applies = hlp == 2;
            bc_dir_1pos_applies = hlp == 1;
            bc_dir_pos_applies = hlp >= 1; 
            bool_v_bc_nodes_applies = sum(this.bool_v_bc_nodes,1) >= 1;
            no_bc = ~bc_dir_pos_applies & ~bool_v_bc_nodes_applies;
            
            if ~isempty(r.DF)
                have_residuals = bc_dir_pos_applies | bool_v_bc_nodes_applies;
                % By sorting and combining the pos/velo Dir BC, the
                % plotting of mixed BCs on one node is plotted correctly.
                residuals_pos = this.bool_u_bc_nodes | this.bool_v_bc_nodes;
                [~, sortidx] = sort([this.idx_u_bc_glob; this.idx_v_bc_glob-posdofs]);
                % Preallocate the residuals matrix
                residuals = zeros(size(residuals_pos));
                maxdfval = max(abs(r.DF(:)))/10;
            end
            
            if r.Vid
                avifile = fullfile(pwd,'output.avi');
                vw = VideoWriter(avifile);
                vw.FrameRate = 25;
                vw.open;
            end
            
            if r.Forces
                % Get forces on each node in x,y,z directions
                % This is where the connection between plane index and
                % x,y,z coordinate is "restored"
                forces = zeros(size(this.bool_u_bc_nodes));
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
                    residual_neumann_forces = zeros(size(this.bool_u_bc_nodes));
                end
            end
            
            %% Loop over time
            ax_extra = this.initRefinedPlot(t,y_dofs,r,pm);
            
            h = pm.nextPlot('geo',sprintf('Deformation at t=%g',t(end)),'x [mm]','y [mm]');
            zlabel(h,'z [mm]');
            
            if ~r.Skel
%                 light('Position',[1 1 1],'Style','infinite','Parent',h);
                musclecol = [0.854688, 0.201563, 0.217188];
            end
            
            axis(h, this.getPlotBox(y_dofs));
            daspect([1 1 1]);
            view(h, [46 30]);
            hold(h,'on');
            
            for ts = 1:length(t)
                % Quit if figure has been closed
                if ~ishandle(h)
                    break;
                end
                u = reshape(y_dofs(1:vstart-1,ts),3,[]);
                v = reshape(y_dofs(vstart:pstart-1,ts),3,[]);
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
                plot3(h,u(1,bool_v_bc_nodes_applies),u(2,bool_v_bc_nodes_applies),u(3,bool_v_bc_nodes_applies),'g.','MarkerSize',20);
                
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
                        forces(1,:),forces(2,:),forces(3,:),0.1,'Color',[.8 .8 1]);
                    % Plot force vectors at face centers
                    for k=1:numfaceswithforce
                        masterfacenodeidx = geo.MasterFaces(force_elem_face_idx(2,k),:);
                        facenodeidx = geo.Elements(force_elem_face_idx(1,k),masterfacenodeidx);
                        
                        facecenter = mean(u(:,facenodeidx),2);
                        quiver3(h,facecenter(1),facecenter(2),facecenter(3),...
                            meanforces(1,k),meanforces(2,k),meanforces(3,k),0.1,'LineWidth',2,'Color','b','MaxHeadSize',1);
                        
                        if ~isempty(r.NF)
                            residual_neumann_forces(this.bc_neum_forces_nodeidx) = r.NF(:,ts);
                            meanforce = mean(residual_neumann_forces(:,facenodeidx),2);
                            quiver3(h,facecenter(1),facecenter(2),facecenter(3),...
                            meanforce(1),meanforce(2),meanforce(3),0.1,'LineWidth',2,'Color','k','MaxHeadSize',1);
                        end
                    end
                end
                
                %% Pressure
                if r.Pressure
                    p = y_dofs(pstart:this.num_uvp_dof,ts);
                    pn = this.idx_p_to_u_nodes;
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
                        u = y_dofs(1:vstart-1,ts);
                        u = u(this.idx_u_glob_elems(:,:,m));
                        
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

                % Call a subroutine for further plots at time step ts
                this.refinedPlot(ax_extra, t, y_dofs, r, ts);
                
                if r.Vid
                    vw.writeVideo(getframe(gcf));
                else
%                     pause(.01);
                      drawnow;
%                     pause;
                end
            end
            
            if r.Vid
                vw.close;
            end

            if isempty(r.PM)
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
            %% Re-add the dirichlet nodes
            uvwall = zeros(this.num_uvp_glob, size(uvw,2));
            uvwall(this.idx_uv_dof_glob,:) = uvw;
            uvwall(this.idx_uv_bc_glob,:) = repmat(this.val_uv_bc_glob,1,size(uvw,2));
            
            sys = this.Model.System;
            zerovel = t > sys.ApplyVelocityBCUntil;
            uvwall(sys.idx_v_bc_glob,zerovel) = 0;
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
        
        function set.UseCrossFibreStiffness(this, value)
            if ~islogical(value) || ~isscalar(value)
                error('UseCrossFibreStiffness must be true or false');
            end
            if this.fUseCrossFibreStiffness ~= value
                this.fUseCrossFibreStiffness = value;
                if ~isempty(this.Model.Config)
                    this.inita0;
                end
            end
        end
        
        function value = get.UseCrossFibreStiffness(this)
            value = this.fUseCrossFibreStiffness;
        end
    end
    
    methods(Access=protected)
        
        function h = initRefinedPlot(this, t, y, r, pm)
            mc = this.Model.Config;
            h = [];
            if ~isempty(mc.Pool) && r.Pool
                h = pm.nextPlot('force','Activation force','t [ms]','alpha');
                axis(h,[0 t(end) 0 1]);
                hold(h,'on');
            end
        end
        
        function refinedPlot(this, h, t, y, r, ts)
            mc = this.Model.Config;
            if ~isempty(mc.Pool) && r.Pool
                dt = this.Model.dt;
                times = 0:dt:ts*dt;
                alpha = mc.Pool.getActivation(times);
                walpha = mc.FibreTypeWeights(1,:,1) * alpha;
                cla(h);
                plot(h,times,alpha);
                plot(h,times,walpha,'LineWidth',2);
            end
        end
        
        function updateDofNums(this, mc)
            tq = mc.PosFE;
            geo_uv = tq.Geometry;
            tl = mc.PressFE;
            geo_p = tl.Geometry;
            this.num_uvp_glob = geo_uv.NumNodes * 6 + geo_p.NumNodes;
            this.num_uvp_dof = this.num_uvp_glob - length(this.val_uv_bc_glob);
        end
        
        function x0 = assembleX0(this)
            % Constant initial values as current node positions
            mc = this.Model.Config;
            tq = mc.PosFE;
            geo = tq.Geometry;
            
            % All zero, especially the first tq.NumNodes * 3 velocity
            % entries and pressure
            x0 = zeros(this.num_uvp_glob,1);
            
            % Fill in the reference configuration positions as initial
            % conditions
            for m = 1:geo.NumElements
                 dofpos = this.idx_u_glob_elems(:,:,m);
                 x0(dofpos) = geo.Nodes(:,geo.Elements(m,:));
            end
            
            % Initial conditions for pressure (s.t. S(X,0) = 0)
            x0(geo.NumNodes*6+1 :end) = -2*this.f.c10-4*this.f.c01;

            % Give the model config a chance to provide x0
            x0 = mc.getX0(x0);

            % Remove dirichlet values
            x0(this.idx_uv_bc_glob) = [];
        end
        
        function MM = assembleMassMatrix(this)
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
            MM(this.idx_uv_bc_glob,:) = [];
            MM(:,this.idx_uv_bc_glob) = [];
            
            % See description of property
            if this.UseDirectMassInversion
                this.Minv = inv(MM(this.idx_v_dof_glob,this.idx_v_dof_glob));
                MM = sparse(I(:),J(:),s(S(:)),3*nd,3*nd);
                % Use identity on left hand side
                MM = blkdiag(speye(size(MM)),...
                    speye(size(MM)),sparse(gp.NumNodes,gp.NumNodes));
                MM(this.idx_uv_bc_glob,:) = [];
                MM(:,this.idx_uv_bc_glob) = [];
            end
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
            D(this.idx_uv_bc_glob,:) = [];
            D(:,this.idx_uv_bc_glob) = [];
            
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
            num_u_glob = geo.NumNodes * 3;
            
            %% Displacement
            this.bool_u_bc_nodes = pos_dir;
            % Set values to node positions
            this.val_u_bc = geo.Nodes(pos_dir);
            this.idx_u_bc_glob = int32(find(pos_dir(:)));
            
            %% Velocity
            % Add any user-defines values (cannot conflict with position
            % dirichlet conditions, this is checked in AModelConfig.getBC)
            this.bool_v_bc_nodes = velo_dir;
            this.idx_v_bc_glob = int32(num_u_glob + find(velo_dir(:)));
            this.val_v_bc = velo_dir_val(velo_dir);
            
            % Compile the global dirichlet values index and value vectors.
            % Here we add zero velocities for each point with fixed position, too.
            [this.idx_uv_bc_glob, sortidx] = sort([this.idx_u_bc_glob; this.idx_u_bc_glob+num_u_glob; this.idx_v_bc_glob]);
            alldirvals = [this.val_u_bc; zeros(size(this.val_u_bc)); this.val_v_bc];
            this.val_uv_bc_glob = alldirvals(sortidx);
            
            % Compute dof positions in global state space vector
            total = geo.NumNodes * 6 + pgeo.NumNodes;
            pos = false(1,total);
            pos(1:num_u_glob) = true;
            pos(this.idx_uv_bc_glob) = [];
            this.idx_u_dof_glob = int32(find(pos));
            this.num_u_dof = length(this.idx_u_dof_glob);
            
            pos = false(1,total);
            pos(num_u_glob+1:num_u_glob*2) = true;
            pos(this.idx_uv_bc_glob) = [];
            this.idx_v_dof_glob = int32(find(pos));
            this.num_v_dof = length(this.idx_v_dof_glob);
            
            % no possible dirichlet values for p. so have all
            this.idx_p_dof_glob = geo.NumNodes*6-length(this.idx_uv_bc_glob) + (1:pgeo.NumNodes);
            this.num_p_dof = length(this.idx_p_dof_glob);
            
            idx = int32(1:total);
            idx(this.idx_uv_bc_glob) = [];
            this.idx_uv_dof_glob = idx;
        end
        
        function [force, nodeidx] = getSpatialExternalForces(this)
            mc = this.Model.Config;
            fe_displ = mc.PosFE;
            geo = fe_displ.Geometry;
            ngp = fe_displ.GaussPointsPerElemFace;
            force = zeros(geo.NumNodes * 3,1);
            faceswithforce = false(1,geo.NumFaces);
            
            globalcoord = strcmp(mc.NeumannCoordinateSystem,'global');
            for fn = 1:geo.NumFaces
                elemidx = geo.Faces(1,fn);
                faceidx = geo.Faces(2,fn);
                masterfacenodeidx = geo.MasterFaces(faceidx,:);
                % So far: Constant pressure on all gauss points!
                P = mc.getBoundaryPressure(elemidx, faceidx);
                if ~isempty(P)
                    faceswithforce(fn) = true;
                    integrand = zeros(3,geo.NodesPerFace);
                    if globalcoord
                        N = geo.FaceNormals(:,faceidx);
                    end
                    for gi = 1:ngp
                        if ~globalcoord
                            N = fe_displ.NormalsOnFaceGP(:,gi,fn);
                        end
                        PN = (P * N) * fe_displ.Ngpface(:,gi,fn)';
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
            ismastercoord = strcmp(mc.a0CoordinateSystem,'master');
            this.HasFibres = false;
            if any(anull(:))
                this.HasFibres = true;
            
                this.a0 = anull;

                % Precomputations
                dNgp = fe.gradN(fe.GaussPoints);
                Ngp = fe.N(fe.GaussPoints);
                
                anulldyadanull = zeros(3,3,fe.GaussPointsPerElem*geo.NumElements);
                dNanull = zeros(geo.DofsPerElement,fe.GaussPointsPerElem,geo.NumElements);
                if this.fUseCrossFibreStiffness
                    a0a0n1 = anulldyadanull;
                    a0a0n2 = anulldyadanull;
                    a0base = anulldyadanull;
                end
                for m = 1 : geo.NumElements
                    u = geo.Nodes(:,geo.Elements(m,:));
                    for gp = 1 : fe.GaussPointsPerElem
                        
                        % Transform a0 to local fibre direction
                        pos = [0 fe.GaussPointsPerElem 2*fe.GaussPointsPerElem]+gp;
                        Jac = u*dNgp(:,pos);
                        if ismastercoord
                            loc_anull = Jac * anull(:,gp,m);
                        else
                            loc_anull = anull(:,gp,m);
                        end
                        loc_anull = loc_anull/norm(loc_anull);
                        
                        % Compute "the" two normals (any will do)
                        loc_anulln1 = circshift(loc_anull,1);
                        loc_anulln1 = loc_anulln1 - (loc_anulln1'*loc_anull) * loc_anull;
                        loc_anulln1 = loc_anulln1/norm(loc_anulln1);

                        loc_anulln2 = circshift(loc_anull,2);
                        loc_anulln2 = loc_anulln2 - (loc_anulln2'*loc_anull) * loc_anull - (loc_anulln2'*loc_anulln1) * loc_anulln1;
                        loc_anulln2 = loc_anulln2/norm(loc_anulln2);
                        
                        % forward transformation of a0 at gauss points
                        % (plotting only so far)
                        if ismastercoord
                            dNanull(:,gp,m) = dNgp(:,pos) * anull(:,gp,m);
                        else
                            dNanull(:,gp,m) = dNgp(:,pos) * (Jac \ anull(:,gp,m));
                        end
                        
                        % a0 dyad a0
                        pos = (m-1)*fe.GaussPointsPerElem+gp;
                        anulldyadanull(:,:,pos) = loc_anull*loc_anull';
                        a0base(:,:,pos) = [loc_anull loc_anulln1 loc_anulln2];
                        if this.fUseCrossFibreStiffness
                            a0a0n1(:,:,pos) = loc_anulln1*loc_anulln1';
                            a0a0n2(:,:,pos) = loc_anulln2*loc_anulln2';
                        end
                    end
                end
                this.a0oa0 = anulldyadanull;
                this.dNa0 = dNanull;
                this.a0Base = a0base;
                
                this.a0oa0n1 = [];
                this.a0oa0n2 = [];
                if this.fUseCrossFibreStiffness
                    this.a0oa0n1 = a0a0n1;
                    this.a0oa0n2 = a0a0n2;    
                end
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
