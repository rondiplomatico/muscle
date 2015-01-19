classdef LongForceBC < muscle.AModelConfig
    % Demo class with a long beam, diagonal fibre direction and two-point
    % boundary face forces in opposing directions.
    
    methods
        function this = LongForceBC
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid([0 20],-40:10:40,[0 15],.1);
            geo = geometry.Cube8Node(pts, cubes);
            this = this@muscle.AModelConfig(geo.toCube20Node);
        end
        
        function configureModel(this, m)
            m.T = 40;
            m.dt = .2;
            m.DefaultMu = [.1; 0; 1; 0; 2.756e-5; 43.373; 7.99; 16.6];
            m.DefaultInput = 1;
            f = m.System.f;
            f.alpha = this.getAlphaRamp(5,.4);
            os = m.ODESolver;
            os.RelTol = .001;
            os.AbsTol = .05;
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
    
end

