classdef Config_FEDEIM < muscle.AModelConfig

    properties(Constant)
        width = 200;
        length = 800;
        height = 200;
    end
    
    methods
        function this = Config_FEDEIM
            elemsize = 50;
            x = -Config_FEDEIM.width/2:elemsize:Config_FEDEIM.width/2;
            y = -Config_FEDEIM.length/2:elemsize:Config_FEDEIM.length/2;
            z = -Config_FEDEIM.height/2:elemsize:Config_FEDEIM.height/2;
            [pts, cubes] = geometry.Cube20Node.DemoGrid(x,y,z);
            geo = geometry.Cube20Node(pts, cubes);
            this = this@muscle.AModelConfig(geo);
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            m.T = 100;
            m.dt = .01;
            m.DefaultMu(1:4) = [1; 20; 0; 0];
            m.DefaultMu(13) = 200;
        end
        
        function prepareSimulation(this, mu, ~)
            sys = this.Model.System;
            sys.f.alpha = this.getAlphaRamp(mu(2),1);
        end
    end
    
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            for k = (0:3)*64+1
                displ_dir(:,geo.Elements(k:k+3,geo.MasterFaces(3,:))) = true;
            end
        end
        
        function anull = seta0(this, anull)
            fe = this.PosFE;
            g = fe.Geometry;
            
            for elem = 1:g.NumElements
                % Get position of Gauss points
                gps = fe.getGlobalGaussPoints(elem);
                fac = gps(2,:)/Config_FEDEIM.length/2;
                maxdeg = pi;
                loc_anull = [-fac .* sin(gps(1,:)/Config_FEDEIM.width*maxdeg);...
                            .5*ones(1,size(gps,2));...
                            -fac .* sin(gps(3,:)/Config_FEDEIM.height*maxdeg)];
                anull(:,:,elem) = loc_anull;
            end
        end
    end
    
end

