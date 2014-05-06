classdef System < models.BaseDynSystem
% MuscleFibreSystem: The global dynamical system used within the MuscleFibreModel
%
% Contains the FibreDynamics, Linear Diffusion and InputConv components as well as a
% ConstInitialValue.
%
% @author Daniel Wirtz @date 2012-11-22
%
% @new{0,7,dw,2012-11-22} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing

    properties
       globidx_displ;
       globidx_pressure;
    end
    
    properties(SetAccess=private)
       DisplFE;
       PressureFE;
       pressure_to_displ_nodes;
       
       % The overall values of dirichlet conditions
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
       bc_dir_displ_val;
       
       bc_dir_velo;
       bc_dir_velo_idx;
       bc_dir_velo_val;
       
       % A helper array containing the indices of the actual degrees of
       % freedom in the global dof indexing.
       %
       % Used to join dofs with dirichlet values
       dof_idx;
    end
    
    methods
        function this = System(model)
            % The system x'' + K(x) is transformed into the first order system
            % v' + K(x), x' = v
            
            % Call superclass constructor
            this = this@models.BaseDynSystem(model);
            
            geo = model.Geometry;
%             tq = triquadratic(geo);
            tq = trilinear(geo);
            tl = trilinear(geo);
    
            % Save for Dynamics
            this.DisplFE = tq;
            this.PressureFE = tl;
            % Find the indices of the pressure nodes in the displacement
            % nodes geometry (used for plotting)
            this.pressure_to_displ_nodes = Utils.findVecInMatrix(tq.nodes,tl.nodes);
            
            %% Dirichlet conditions: Position (fix one side)
            displ_dir = false(3,tq.NumNodes);
            %% Normal config
            % Fix the front of the first four cubes
%             for k = [5 6 11 12]
            for k = 1
                % Quadratic
%                 displ_dir(:,tq.elems(k,[6:8 11 12 18:20])) = true;
                % Linear
%                displ_dir(:,tq.elems(k,5:8)) = true;
            end

            %% Dirichlet conditions: Velocity (fix one side)
            % Compile indices of velocity DOFS that will be increased over
            % time as BC
            velo_dir = false(3,tq.NumNodes);
            velo_dir_val = zeros(3,tq.NumNodes);
%             for k = [1 2 7 8]
            for k = 1
                % Only z direction
                % Quadratic
%                 velo_dir(3,tq.elems(k,[1:3 9 10 13:15])) = true;
                % Linear
%                 velo_dir(3,tq.elems(k,1:4)) = true;
            end
%             velo_dir_val(velo_dir) = .005;
            
            % Call subroutine for boundary condition index crunching
            this.computeBC(displ_dir, velo_dir, velo_dir_val);
            
            % Construct global indices in uvw from element nodes. Each dof in
            % an element is used three times for x,y,z displacement. The
            % "elems" matrix contains the overall DOF numbers of each
            % element in the order of the nodes (along row) in the master
            % element.
            ne = tq.NumElems;
            globalelementdofs = zeros(3,tq.DofsPerElement,ne);
            for m = 1:ne
                % First index of element dof in global array
                hlp = (tq.elems(m,:)-1)*3+1;
                % Use first, second and third as positions.
                globalelementdofs(:,:,m) = [hlp; hlp+1; hlp+2];
            end
            this.globidx_displ = globalelementdofs;
            
            % The same for the pressure
            globalpressuredofs = zeros(tl.DofsPerElement,tl.NumElems);
            for m = 1:tl.NumElems
                globalpressuredofs(:,m) = tq.NumNodes * 6 + tl.elems(m,:);
            end
            this.globidx_pressure = globalpressuredofs;
            
            %% Compile Mass Matrix
            % Augment mass matrix for all 3 displacement directions
            nd = tq.NumNodes;
            [i, j, s] = find(tq.M);
            I = [3*(i'-1)+1; 3*(i'-1)+2; 3*(i'-1)+3];
            J = [3*(j'-1)+1; 3*(j'-1)+2; 3*(j'-1)+3];
            S = repmat(1:length(s),3,1);
            MM = sparse(I(:),J(:),s(S(:)),3*nd,3*nd);
            % Insert identity for velocity and zeros for pressure
            MM = blkdiag(speye(size(MM)),...
                MM,sparse(tl.NumNodes,tl.NumNodes));
            
            % Strip out the entries of dirichlet nodes
            MM(this.bc_dir_idx,:) = [];
            MM(:,this.bc_dir_idx) = [];
            
            this.M = dscomponents.ConstMassMatrix(MM);
            
            %% Set system components
            % Core nonlinearity
            this.f = muscle.Dynamics(this);
            
            %% Initial value
            this.x0 = this.assembleX0;
        end
        
        function pm = plot(this, t, uvw, pm)
            if nargin < 4
                pm = PlotManager;
                pm.LeaveOpen = true;
            end
            
            dfem = this.DisplFE;
            pfem = this.PressureFE;
            vstart = dfem.NumNodes * 3+1;
            pstart = dfem.NumNodes * 6+1;
            e = dfem.edges;
            
            %% Re-add the dirichlet nodes
            yall = zeros(dfem.NumNodes * 6 + pfem.NumNodes, size(uvw,2));
            yall(this.dof_idx,:) = uvw;
            for k=1:length(this.bc_dir_idx)
                yall(this.bc_dir_idx(k),:) = this.bc_dir_val(k);
            end
            uvw = yall;
            
            bc_dir_displ_applies = sum(this.bc_dir_displ,1) >= 1;
            bc_dir_velo_applies = sum(this.bc_dir_velo,1) >= 1;
            
            %% Loop over time
            h = pm.nextPlot('geo','Output','x','y');
            for ts = 1:1:length(t)
                % Quit if figure has been closed
                if ~ishandle(h)
                    break;
                end
                u = reshape(uvw(1:vstart-1,ts),3,[]);
                v = reshape(uvw(vstart:pstart-1,ts),3,[]);
                
                plot3(h,u(1,:),u(2,:),u(3,:),'k.','MarkerSize',14);
                hold(h,'on');
                
                % Velocities
                %quiver3(h,u(1,:),u(2,:),u(3,:),v(1,:),v(2,:),v(3,:),'k.','MarkerSize',14);
                quiver3(h,u(1,:),u(2,:),u(3,:),v(1,:),v(2,:),v(3,:),'MarkerSize',14);
                
                %% Dirichlet conditions
                % Displacement
                plot3(h,u(1,bc_dir_displ_applies),u(2,bc_dir_displ_applies),u(3,bc_dir_displ_applies),'k.','MarkerSize',20);
                for k=1:size(e,1)
                    plot3(h,u(1,[e(k,1) e(k,2)]),u(2,[e(k,1) e(k,2)]),u(3,[e(k,1) e(k,2)]),'r');
                end
                
                % Velocity
                plot3(h,u(1,bc_dir_velo_applies),u(2,bc_dir_velo_applies),u(3,bc_dir_velo_applies),'g.','MarkerSize',20);
                quiver3(h,u(1,:),u(2,:),u(3,:),v(1,:),v(2,:),v(3,:),'MarkerSize',14);
%                 for k=1:size(e,1)
%                     plot3(h,u(1,[e(k,1) e(k,2)]),u(2,[e(k,1) e(k,2)]),u(3,[e(k,1) e(k,2)]),'r');
%                 end
                
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
                
                %% Misc
                axis(h,'tight');
                view(h, [46 30]);
                title(h,sprintf('Deformation at t=%g',t(ts)));
                hold(h,'off');
                
%                 pause(1);
                pause;
            end

            if nargin < 4
                pm.done;
            end
        end
    end
    
    methods(Access=private)
        function x0 = assembleX0(this)
            % Constant initial values as current node positions
            tq = this.DisplFE;
            tl = this.PressureFE;
            % All zero, especially the first tq.NumNodes * 3 velocity
            % entries and pressure
            x0 = zeros(tq.NumNodes * 6 + tl.NumNodes,1);
            
            % Fill in the reference configuration positions as initial
            % conditions
            for m = 1:tq.NumElems
                 dofpos = this.globidx_displ(:,:,m);
                 x0(dofpos) = tq.nodes(:,tq.elems(m,:));
            end
            
            % Initial conditions for pressure (s.t. S(X,0) = 0)
            x0(tq.NumNodes * 6+1:end) = -2*this.f.c10-4*this.f.c01;
            
            % Remove dirichlet values
            x0(this.bc_dir_idx) = [];
            
            x0 = dscomponents.ConstInitialValue(x0);
        end
        
        function computeBC(this, displ_dir, velo_dir, velo_dir_val)
            
            if any(any(displ_dir & velo_dir))
                error('Cannot impose displacement and velocity dirichlet conditions on same DoF');
            end
            fe_displ = this.DisplFE;
            fe_press = this.PressureFE;
            
            % Position of position entries in global state space vector
            velocitydofs_start_global = fe_displ.NumNodes * 3;
            
            %% Displacement
            this.bc_dir_displ = displ_dir;
            % Set values to node positions
            this.bc_dir_displ_val = fe_displ.nodes(displ_dir);
            % Add zeros for respective velocities
            this.bc_dir_displ_val = [this.bc_dir_displ_val; zeros(size(this.bc_dir_displ_val))];
            % Collect indices of dirichlet values per 3-set of x,y,z values
            relpos = find(displ_dir(:));
            % Same positions for points and velocity
            this.bc_dir_displ_idx = [relpos; relpos + velocitydofs_start_global];
            
            %% Velocity
            this.bc_dir_velo = velo_dir;
            this.bc_dir_velo_val = velo_dir_val(velo_dir);
            this.bc_dir_velo_idx = velocitydofs_start_global + find(velo_dir(:));
            
            %% Hydrostatic Pressure Dirichlet conditions
%             press_dir = false(1,tl.NumNodes);
%             press_dir(1) = true;
%             this.bc_dir = [displ_dir; press_dir];

            % The according velocities are all zero for dirichlet points
            this.bc_dir_val = [this.bc_dir_displ_val; this.bc_dir_velo_val];
            this.bc_dir_idx = [this.bc_dir_displ_idx; this.bc_dir_velo_idx];
            
            % Compute dof positions in global state space vector
            pos = 1:(fe_displ.NumNodes * 6 + fe_press.NumNodes);
            pos(this.bc_dir_idx) = [];
            this.dof_idx = pos;
        end
    end
    
end
