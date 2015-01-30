classdef Cube12 < muscle.AModelConfig
    
    methods
        function this = Cube12
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid(-1:1,-1:2,-1:1);
            
            geo = geometry.Cube8Node(pts, cubes);
            this = this@muscle.AModelConfig(geo.toCube20Node);
        end
        
        function configureModel(this, m)
            m.T = 40;
            m.dt = .05;
            m.DefaultMu(2) = 20;
            model.DefaultMu(13) = 200; % [kPa]
        end
    end
    
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            for k = [5 6 11 12]
                displ_dir(:,geo.Elements(k,geo.MasterFaces(4,:))) = true;
            end
        end
        
%         function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
%             %% Dirichlet conditions: Position (fix one side)
%             geo = this.PosFE.Geometry;
%             for k = [1 2 7 8]  
%                 velo_dir(1,geo.Elements(k,geo.MasterFaces(3,:))) = true;
%             end
%             velo_dir_val(velo_dir) = -.1;
%         end
        
        function anull = seta0(~, anull)
            % Direction is xz
            anull([1 3],:,:) = 1;
        end
    end
    
end

