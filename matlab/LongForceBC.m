classdef LongForceBC < muscle.AModelConfig
    
    methods
        function this = LongForceBC
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid([0 20],-10:10:10,[0 15],.1);
            geo = geometry.Cube8Node(pts, cubes);
            this = this@muscle.AModelConfig(geo);
        end
        
        function configureModel(~, model)
            model.T = 5;
            model.dt = .1;
            f = model.System.f;
            f.alpha = .3;
            f.Viscosity = 0;
            os = model.ODESolver;
            os.RelTol = .1;
            os.AbsTol = .1;
            model.System.Inputs{1} = @(t)min(1,t);
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
           anull([1 3],:,:) = 1;
        end
    end
    
end

