classdef Long < muscle.AModelConfig
% A long geometry with 20% deviation from default cubic positions and
% complex fibre structure

    methods
        function this = Long(devi)
            if nargin < 1
                devi = .2;
            end
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube20Node.DemoGrid(-10:10:10,-40:10:40, [0 10], devi);
            geo = geometry.Cube20Node(pts, cubes);
            this = this@muscle.AModelConfig(geo);
        end
        
        function configureModel(this, m)
            m.T = 50;
            m.dt = 1;
            m.DefaultMu = [.1; 0; 0; 0];
            f = m.System.f;
            f.alpha = this.getAlphaRamp(10,1);
        end
    end
    
    methods(Access=protected)
        
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
    
end

