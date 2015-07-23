classdef CubePull < muscle.AModelConfig
    
    methods
        function this = CubePull(varargin)
            % Single cube with same config as reference element
            this = this@muscle.AModelConfig(varargin{:});
            this.init;
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            m.T = 150;
            m.dt = .5;
            m.DefaultMu(1) = .5;
            m.DefaultInput = 1;
        end
        
        function configureModelFinal(this)
            % Set desired reduction to full state space dimension by
            % default
            m = this.Model;
            m.SpaceReducer.Value = m.System.num_u_dof;
            configureModelFinal@muscle.AModelConfig(this);
        end
        
        function P = getBoundaryPressure(~, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is megaPascal [MPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            if elemidx == 1 && faceidx == 3
                P = 1;
            end
        end
        
        function u = getInputs(this)
            u = {this.getAlphaRamp(this.Model.T/3,1)};
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(~)
            %[pts, cubes] = geometry.Cube8Node.DemoGrid([0 2.5 5],-1:20,[0 2.5 5]);
            [pts, cubes] = geometry.Cube8Node.DemoGrid([0 5],-1:10,[0 5]);
            geo = geometry.Cube8Node(pts, cubes);
            geo = geo.toCube27Node;
            %geo = geo.toCube20Node;
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            % Fix all on left and only the y,z directions of the back right
            % nodes
            displ_dir(2,geo.Elements(end,geo.MasterFaces(4,:))) = true;
            displ_dir(:,geo.Elements(end,geo.MasterFaces(4,1))) = true;
        end
        
        function anull = seta0(~, anull)
           % Direction is xz
           %anull([1 3],:,:) = 1;
        end
    end
    
    methods(Static)
        function test_CubePull
            m = muscle.Model(CubePull);
            mu = m.getRandomParam;
            mu(3) = .1;
            m.simulateAndPlot(true,mu,1);
        end
    end
end

