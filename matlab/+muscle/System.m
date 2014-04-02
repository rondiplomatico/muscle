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
    end
    
    properties(Transient)
       
    end
    
    properties(SetAccess=private)
       
    end
    
    properties (Constant)        
       
    end
    
    properties(Access=private)
       
    end
    
    methods
        function this = System(model)
            % Call superclass constructor
            this = this@models.BaseDynSystem(model);
            
            geo = model.Geometry;
            %tq = triquadratic(geo);
            tq = trilinear(geo);
            % Save for Dynamics
            this.DisplFE = tq;
            
            % Construct global indices in y from element dofs. Each dof in
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
            this.globidx = globalelementdofs;
            
            % Linear diffusion part
%             this.A = this.assembleA;

            %% Get Mass Matrix
            % Augment mass matrix for all 3 displacement directions
            nd = tq.NumDofs;
            [i, j, s] = find(tq.M);
            I = [3*(i'-1)+1; 3*(i'-1)+2; 3*(i'-1)+3];
            J = [3*(j'-1)+1; 3*(j'-1)+2; 3*(j'-1)+3];
            S = repmat(1:length(s),3,1);
            M3 = sparse(I(:),J(:),s(S(:)),3*nd,3*nd);
            
            this.M = dscomponents.ConstMassMatrix(M3);
            %this.M = dscomponents.ConstMassMatrix(blkdiag(speye(size(M3)),M3));
            
            %% Set system components
            % Core nonlinearity
            this.f = muscle.Dynamics(this);
            
%             % Linear input B for motoneuron
%             this.B = this.assembleB;
%             
            this.x0 = this.assembleX0;
        end
        
        function pm = plot(this, t, y, pm)
            if nargin < 4
                pm = PlotManager(false,2,2);
                pm.LeaveOpen = true;
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
            x0 = zeros(tq.NumDofs * 3,1);
            for m = 1:tq.NumElems
                 dofpos = this.globidx(:,:,m);
                 x0(dofpos) = tq.dofs(:,tq.elems(m,:));
            end
            x0 = dscomponents.ConstInitialValue(x0);
        end
    end
    
end
