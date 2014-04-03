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
       DisplFE;
       
       globidx;
       
       % Boundary conditions: 3 times numnodes logical matrix telling
       % which degree of freedom of which (position) node is fixed
       %     n1   n2   n3 ...
       % x | 1    0    0       |
       % y | 1    1    0 ...   |
       % z | 1    0    0       |
       bc_dir;
       
       % The dirichlet values (linear array)
       bc_dir_val;
       
       % The positions of the dirichlet values within the GLOBAL node/xyz
       % state space vector
       bc_dir_idx;
    end
    
    methods
        function this = System(model)
            % The system x'' + K(x) is transformed into the first order system
            % v' + K(x), x' = v
            
            % Call superclass constructor
            this = this@models.BaseDynSystem(model);
            
            geo = model.Geometry;
            tq = triquadratic(geo);
            %tq = trilinear(geo);
            % Save for Dynamics
            this.DisplFE = tq;
            
            %% Dirichlet conditions
            dir = false(3,tq.NumNodes);
            dir(:,1) = true;
            dir(:,2) = true;
            %this.bc_dir = dir;
            
            % Position of position entries in global state space vector
            positiondofs_start_global = tq.NumNodes * 3;
            
            % The according velocities are all zero for dirichlet points
            this.bc_dir_val = [zeros(size(tq.nodes(dir))); tq.nodes(dir)];
            
            % Collect indices of dirichlet values per 3-set of x,y,z values
            bc_dir_relpos = find(dir(:));
            % Same positions for points and velocity
            this.bc_dir_idx = [bc_dir_relpos; bc_dir_relpos + positiondofs_start_global];
            
            % Construct global indices in y from element nodes. Each dof in
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
            % The element dofs are for the positions and not velocities
            globalelementdofs = globalelementdofs + positiondofs_start_global;
            this.globidx = globalelementdofs;
            
            %% Compile Mass Matrix
            % Augment mass matrix for all 3 displacement directions
            nd = tq.NumNodes;
            [i, j, s] = find(tq.M);
            I = [3*(i'-1)+1; 3*(i'-1)+2; 3*(i'-1)+3];
            J = [3*(j'-1)+1; 3*(j'-1)+2; 3*(j'-1)+3];
            S = repmat(1:length(s),3,1);
            M3 = sparse(I(:),J(:),s(S(:)),3*nd,3*nd);
            
            % Strip out the entries of dirichlet nodes:
            % This is done before adding the eye-copy for the velocities,
            % as they are reduced by the same amount of dofs.
            M3(bc_dir_relpos,:) = [];
            M3(:,bc_dir_relpos) = [];
            
            % Hence, the mass matrix is for the v' part, whereas the x'
            % part has simple identity
            this.M = dscomponents.ConstMassMatrix(blkdiag(speye(size(M3)),M3));
            
            %% Set system components
            % Core nonlinearity
            this.f = muscle.Dynamics(this);
            
            %% Initial value
            this.x0 = this.assembleX0;
        end
        
        function pm = plot(this, t, y, pm)
            if nargin < 4
                pm = PlotManager;
                pm.LeaveOpen = true;
            end
            
            dfem = this.DisplFE;
            xstart = dfem.NumNodes * 3 + 1;
            e = dfem.edges;
            
            h = pm.nextPlot('geo','Output','x','y');
            for ts = 1:length(t)
                
                p = reshape(y(xstart:end,ts),3,[]);
                plot3(h,p(1,:),p(2,:),p(3,:),'k.','MarkerSize',14);
                hold(h,'on');
                for k=1:size(e,1)
                    plot3(h,p(1,[e(k,1) e(k,2)]),p(2,[e(k,1) e(k,2)]),p(3,[e(k,1) e(k,2)]),'r');
                end
                axis(h,'tight');
                title(h,sprintf('Deformation at t=%g',t(ts)));
                hold(h,'off');
                
                pause(.4);
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
            % All zero, especially the first tq.NumNodes * 3 velocity
            % entries
            x0 = zeros(tq.NumNodes * 6,1);
            
            % Fill in the reference configuration positions as initial
            % conditions
            for m = 1:tq.NumElems
                 dofpos = this.globidx(:,:,m);
                 x0(dofpos) = tq.nodes(:,tq.elems(m,:));
            end
            
            % Remove dirichlet values
            x0(this.bc_dir_idx) = [];
            
            x0 = dscomponents.ConstInitialValue(x0);
        end
    end
    
end
