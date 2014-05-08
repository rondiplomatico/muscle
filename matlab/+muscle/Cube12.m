classdef Cube12 < muscle.AModelConfig
    
    methods
        function this = Cube12
            % Single cube with same config as reference element
            [pts, cubes] = cubegeom.DemoCubeGrid(-1:1,-1:2,-1:1);
            
            geo = cubegeom(pts, cubes);
            this = this@muscle.AModelConfig(geo);
        end
    end
    
    methods(Access=protected)
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            tq = this.PosFE;
            for k = [5 6 11 12]
                displ_dir(:,tq.elems(k,[6:8 11 12 18:20])) = true;
            end
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            %% Dirichlet conditions: Position (fix one side)
            tq = this.PosFE;
            
            for k = [1 2 7 8]  
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
        
        function anull = seta0(~, anull)
            % Direction is x
            anull(2,:,:) = -1;
            anull(1,:,:) = 1;
        end
    end
    
end

