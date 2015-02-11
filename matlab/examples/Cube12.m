classdef Cube12 < muscle.AModelConfig
    
    methods
        function this = Cube12(varargin)
            this = this@muscle.AModelConfig(varargin{:});
            this.init;
        end
        
        function configureModel(~, m)
            m.T = 40;
            m.dt = .05;
            % Activate over 20ms
            m.DefaultMu(2) = 20;
            m.DefaultMu(13) = 100; % [kPa]
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(~)
            [pts, cubes] = geometry.Cube8Node.DemoGrid(-1:1,-1:2,-1:1);
            geo = geometry.Cube8Node(pts, cubes);
            geo = geo.toCube20Node;
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.PosFE.Geometry;
            for k = [5 6 11 12]
                displ_dir(:,geo.Elements(k,geo.MasterFaces(4,:))) = true;
            end
        end

        function anull = seta0(~, anull)
            % Direction is xz
            anull([1 3],:,:) = 1;
        end
    end
    
    methods(Static)
        function test_Cube12
            m = muscle.Model(Cube12);
            m.simulateAndPlot;
        end
    end
    
end

