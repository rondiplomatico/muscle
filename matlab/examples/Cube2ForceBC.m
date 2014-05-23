classdef Cube2ForceBC < muscle.AModelConfig
    
    methods
        function this = Cube2ForceBC
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid(0:1,-1:2,0:1);
            geo = geometry.Cube8Node(pts, cubes);
            this = this@muscle.AModelConfig(geo.toCube20Node, geo);
        end
        
        function configureModel(~, model)
            model.T = 4;
            model.dt = .01;
            f = model.System.f;
            f.alpha = 0;
            f.Viscosity = .1;
            os = model.ODESolver;
            os.RelTol = .0001;
            os.AbsTol = .05;
            % Ramp up the external pressure
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
%             if elemidx == 1 && faceidx == 3
%                 P = zeros(3);
%                 P(3,2) = 1;
%                 %P(2,2) = 1;
%                 P(1,2) = .5;
%             end
            if elemidx == 1 && faceidx == 2
                P = -.8;
            end
            if elemidx == 1 && faceidx == 5
                P = -.4;
            end
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
%            anull([1 3],:,:) = 1;
        end
    end
    
end

