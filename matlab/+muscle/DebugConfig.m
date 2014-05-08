classdef DebugConfig < muscle.AModelConfig
    
    methods
        function this = DebugConfig
            % Single cube with same config as reference element
            [pts, cubes] = cubegeom.DemoCubeGrid(0:1,0:1,0:1);
            pts(pts == 0) = -1;
            geo = cubegeom(pts, cubes);
            this = this@muscle.AModelConfig(geo);
        end
        
    end
    
    methods(Access=protected)
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            tq = this.PosFE;
            displ_dir(:,tq.elems(1,[6:8 11 12 18:20])) = true;
        end
        
        function anull = seta0(this, anull)
            % Set discrete a0 values at all gauss points
            
            % Direction is x
            anull(2,:,:) = -1;
            anull(1,:,:) = 1;
        end
    end
    
end

