classdef Long < muscle.AModelConfig
% A long geometry with 20% deviation from default cubic positions and
% complex fibre structure

    methods
        function this = Long(varargin)
            this = this@muscle.AModelConfig(varargin{:});
            this.addOption('Devi',1);
            this.init;
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            m.T = 50;
            m.dt = 1;
            m.DefaultMu(1) = .1;
            m.DefaultMu(2) = 20;
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube20Node.DemoGrid(-10:10:10,-40:10:40, [0 10], this.Options.Devi);
            geo = geometry.Cube20Node(pts, cubes);
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            for k = [8 16]
                displ_dir(:,geo.Elements(k,[6:8 11 12 18:20])) = true;
            end
        end
        
        function anull = seta0(this, anull)
            fe = this.PosFE;
            if fe.GaussPointsPerElem ~= 27
                warning('a0 designed for 27 gauss points!');
            end
            x = linspace(0,1,8*3);
            basea0 = [sin(x*pi); cos(x*pi)]; %; zeros(size(x))
            front = fe.GaussPoints(2,:) < 0;
            mid = fe.GaussPoints(2,:) == 0;
            back = fe.GaussPoints(2,:) > 0;
            % Direction is x
            for m = 1:8
                off = (m-1)*3;
                anull([1 3],front,[m m+8]) = basea0(2,off+1);
                anull([1 3],mid,[m m+8]) = basea0(2,off+2);
                anull([1 3],back,[m m+8]) = basea0(2,off+3);
                
                anull(2,front,[m m+8]) = basea0(1,off+1);
                anull(2,mid,[m m+8]) = basea0(1,off+2);
                anull(2,back,[m m+8]) = basea0(1,off+3);
            end
        end
    end
    
    methods(Static)
        function test_Long
            m = muscle.Model(Long);
            m.simulateAndPlot;
        end
    end
    
end

