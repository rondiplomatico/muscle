classdef LongForceBC < muscle.AModelConfig
    % Demo class with a long beam, diagonal fibre direction and two-point
    % boundary face forces in opposing directions.
    
    methods
        function this = LongForceBC(varargin)
            this = this@muscle.AModelConfig(varargin{:});
            this.init;

            this.ActivationRampMax = .4;
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            m.T = 40;
            m.dt = .2;
            m.DefaultMu(1) = .1;
            m.DefaultMu(2) = 5;
            m.DefaultMu(3) = 100;
            m.DefaultInput = 1;
            os = m.ODESolver;
            os.RelTol = .001;
            os.AbsTol = .05;
        end
        
        function u = getInputs(this)
            % Returns the inputs `u(t)` of the model.
            %
            % if neumann boundary conditions are used, this input is
            % multiplied with the mu(3) parameter, which determines the
            % maximum force that is applied. u(t) determines its temporal
            % strength.
            %
            % this.Model can be used to get access to the model this
            % configuration is applied to.
            %
            % Return values:
            % u: The cell array of input functions to use within this
            % model.
            %
            % @type cell @default {}
            u = {this.getAlphaRamp(30,1)};
        end
        
        function P = getBoundaryPressure(~, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is kiloPascal [kPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            if elemidx == 1 && faceidx == 1
                P = [1 0 0
                     0 0 0;
                     0 0 0];
            elseif elemidx == 5 && faceidx == 2
                P = [1 0 0
                     0 0 0;
                    .5 0 0];
            end
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid([0 20],-40:10:40,[0 15],.1);
            geo = geometry.Cube8Node(pts, cubes);
            geo = geo.toCube20Node;
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            displ_dir(:,geo.Elements(8,geo.MasterFaces(4,:))) = true;
        end
        
        function anull = seta0(~, anull)
           % Direction is xz
           anull([1 3],:,:) = 1;
        end
    end
    
    methods(Static)
        function test_LongForceBC
            m = muscle.Model(LongForceBC);
            m.simulateAndPlot;
        end
    end
    
end

