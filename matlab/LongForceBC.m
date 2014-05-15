classdef LongForceBC < muscle.AModelConfig
    
    methods
        function this = LongForceBC
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoCubeGrid(0:1,-1:1,0:1);
            geo = geometry.Cube8Node(pts, cubes);
            this = this@muscle.AModelConfig(geo);
        end
        
        function configureModel(~, model)
            model.T = 50;
            model.dt = 1;
            f = model.System.f;
            f.alpha = 0;
            f.Viscosity = 0;
            os = model.ODESolver;
            os.RelTol = .3;
            os.AbsTol = .1;
        end
        
        function P = getBoundaryPressure(~, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is kiloPascal [kPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            if elemidx == 1 && faceidx == 3
                P = -.1;
            end
        end
    end
    
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            displ_dir(:,geo.Elements(2,[6:8 11 12 18:20])) = true;
        end
        
        function anull = seta0(~, anull)
           % Direction is xz
%            anull([1 3],:,:) = 1;
        end
    end
    
end

