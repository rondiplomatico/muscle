classdef FibreDirections < muscle.AModelConfig
    % An example illustrating the fibre direction options.
    % 
    
    methods
        function this = FibreDirections(version)
            if nargin < 1
                version = 1;
            end
            this = this@muscle.AModelConfig();
            this.init;
            
            if version == 2
                this.a0CoordinateSystem = 'reference';
            else
                this.a0CoordinateSystem = 'master';
            end
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            m.T = 50;
            m.dt = 1;
            m.ODESolver.RelTol = 1e-6;
            m.ODESolver.AbsTol = 1e-6;
            m.DefaultMu(2) = 30;
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(~)
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid([0 1],[0 1],[0 1]);
            theta = -.3;
            R = [cos(theta) -sin(theta) 0
                 sin(theta) cos(theta)  0
                 0          0           1];
            pts = circshift(R,[1 1])*R*pts;
            geo = geometry.Cube8Node(pts, cubes);
            geo = geo.toCube27Node;
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.PosFE.Geometry;
            displ_dir(:,geo.Elements(1,geo.MasterFaces(1,:))) = true;
            %displ_dir(1,geo.Elements(1,geo.MasterFaces(1,:))) = true;
            %displ_dir(:,geo.Elements(1,geo.MasterFaces(1,5))) = true;
        end
        
        function anull = seta0(~, anull)
            anull(1,:,:) = 1;
        end
    end
    
    methods(Static)
        function test_FibreDirections
            f = FibreDirections(1);
            m=f.createModel;
            m.plotGeometrySetup;
            m.simulateAndPlot;
            
            f = FibreDirections(2);
            m=f.createModel;
            m.plotGeometrySetup;
            m.simulateAndPlot;
        end
    end
    
end

