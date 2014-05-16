classdef DebugConfig < muscle.AModelConfig
    % A simple configuration for Debug purposes.
    % 
    % Uses a single undeformed 
    
    properties(SetAccess=private)
        Version;
    end
    
    methods
        function this = DebugConfig(version)
            % Creates a Debug simple muscle model configuration.
            %
            % Version 1: No fibres
            % Version 2: Fibres, but no FibreTypes
            % Version 3: Fibres with fibre type distribution and APExp
            % evaluations
            if nargin < 1
                version = 3;
            end
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid(0:1,0:1,0:1);
            pts(pts == 0) = -1;
            geo = geometry.Cube8Node(pts, cubes);
            this = this@muscle.AModelConfig(geo);
            
            this.Version = version;
            %% Muscle fibre weights
            if version == 3
                types = [0 .2 .4 .6 .8 1];
                ftw = zeros(geo.GaussPointsPerElem,length(types),geo.NumElements);
                % Test: Use only slow-twitch muscles
                ftw(:,1,:) = 1;
                this.FibreTypeWeights = ftw;
                this.FibreTypes = types;
            end
        end
        
    end
    
    methods(Access=protected)
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            %displ_dir(:,geo.Elements(1,[6:8 11 12 18:20])) = true;
            
            % Fix left nodes completely
            displ_dir(:,geo.Elements(1,[6 7 11 18 19])) = true;
            
            % Only fix the right back corner nodes in yz-direction
            displ_dir([2 3],geo.Elements(1,[8 12 20])) = true;
        end
        
        function anull = seta0(this, anull)
            % Set discrete a0 values at all gauss points
            
            if this.Version >= 2
%                 anull(2,:,:) = -1;
                anull(1,:,:) = 1;
            end
        end
    end
    
end

