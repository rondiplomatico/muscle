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
            
            % Linear diffusion part
%             this.A = this.assembleA;

            % Get Mass Matrix
            this.M = dscomponents.ConstMassMatrix(model.Geometry.M);
            
            %% Set system components
            % Core nonlinearity
            this.f = muscle.Dynamics(this);
            
            % Linear input B for motoneuron
            this.B = this.assembleB;
            
            % Constant initial values
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
       
        function M = assembleM(this)
            % input conversion matrix, depends on fibre type. Has only one
            % entry in second row.
            g = this.Model.Geometry;
            done = [];
            for k=1:g.NumPoints
                n = g.neighbors{k};
                
                
            end
            M = dscomponents.ConstMassMatrix(M);
        end
        
    end
    
end
