classdef QuasiStaticTest < muscle.AModelConfig
% A long geometry with 20% deviation from default cubic positions and
% complex fibre structure

    methods
        function this = QuasiStaticTest
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid(0:20:40,[0 20],[0 20]);
            geo = geometry.Cube8Node(pts, cubes);
            this = this@muscle.AModelConfig(geo);
        end
        
        function configureModel(~, m)
            
        end
    end
    
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            % Fix front and back
            displ_dir(:,geo.Elements(1,geo.MasterFaces(1,:))) = true;
        end
        
        function anull = seta0(this, anull)
            % Fibres in x direction
            anull(1,:,:) = 1;
        end
    end
    
end

