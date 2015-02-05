classdef Cube2ForceBC < muscle.AModelConfig
    
    methods
        function this = Cube2ForceBC
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid(0:1,-1:2,0:1);
            geo = geometry.Cube8Node(pts, cubes);
            this = this@muscle.AModelConfig(geo.toCube20Node);
        end
        
        function configureModel(this, m)
            m.T = 4;
            m.dt = .01;
            m.DefaultMu(1) = .1;
            m.DefaultMu(3) = 100;
            m.DefaultInput = 1;
            os = m.ODESolver;
            os.RelTol = .0001;
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
            if elemidx == 1 && faceidx == 2
                P = -.8;
            end
            if elemidx == 1 && faceidx == 5
                P = -.4;
            end
        end
        
        function u = getInputs(this)
            u = {this.getAlphaRamp(1,1)};
        end
    end
    
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            % Fix all on left and only the y,z directions of the back right
            % nodes
            displ_dir(:,geo.Elements(3,geo.MasterFaces(4,:))) = true;
            %displ_dir(1,geo.Elements(3,[9 18 27])) = false;
            displ_dir(1,geo.Elements(3,[8 16 20])) = false;
        end
        
        function anull = seta0(~, anull)
           % Direction is xz
           anull([1 3],:,:) = 1;
        end
    end
    
end

