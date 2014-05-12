classdef Long < muscle.AModelConfig
    
    methods
        function this = Long
            % Single cube with same config as reference element
            [pts, cubes] = cubegeom.DemoCubeGrid(-1:1,-4:4,0:1);
            
            geo = cubegeom(pts, cubes);
            this = this@muscle.AModelConfig(geo);
        end
    end
    
    methods(Access=protected)
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            tq = this.PosFE;
            for k = [8 16]
                displ_dir(:,tq.elems(k,[6:8 11 12 18:20])) = true;
            end
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            %% Dirichlet conditions: Position (fix one side)
            tq = this.PosFE;
            
            for k = [1 9]
                % Quadratic
                velo_dir(1,tq.elems(k,[1:3 9 10 13:15])) = true;
%                 velo_dir(1,tq.elems(k,[3 10 15])) = true;
            end
            velo_dir_val(velo_dir) = -.1;

%             velo_dir(1,tq.elems(1,[1 2 9 13 14])) = true;
%             velo_dir(1,tq.elems(2,[1 2 9 13 14])) = true;
%             velo_dir(1,tq.elems(7,[2 3 10 14 15])) = true;
%             velo_dir(1,tq.elems(8,[2 3 10 14 15])) = true;
%             velo_dir_val(1,tq.elems(1,[1 2 9 13 14])) = -.01;
%             velo_dir_val(1,tq.elems(2,[1 2 9 13 14])) = -.01;
%             velo_dir_val(2,tq.elems(7,[2 3 10 14 15])) = .01;
%             velo_dir_val(2,tq.elems(8,[2 3 10 14 15])) = .01;
        end
        
        function anull = seta0(this, anull)
            geo = this.Geometry;
            if geo.NumGaussp ~= 27
                warning('a0 designed for 27 gauss points!');
            end
            x = linspace(0,1,8*3);
            basea0 = [sin(x*pi); cos(x*pi)]; %; zeros(size(x))
            front = geo.gaussp(2,:) < 0;
            mid = geo.gaussp(2,:) == 0;
            back = geo.gaussp(2,:) > 0;
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

