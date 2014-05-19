classdef DebugConfig < muscle.AModelConfig
    % A simple configuration for Debug purposes.
    % 
    % Uses a single undeformed cube with triquadratic position shape
    % functions and trilinear shape functions for pressure.
    %
    %
    % Version 1: No fibres, and one fixed side
    % Version 2: Fibres, but no FibreTypes
    % Version 3: Fibres with fibre type distribution and APExp
    % evaluations
    
    properties(SetAccess=private)
        Version;
    end
    
    methods
        function this = DebugConfig(version)
            % Creates a Debug simple muscle model configuration.
            %
            
            if nargin < 1
                version = 3;
            end
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube8Node.DemoGrid([-1 1],[-1 1],[-1 1]);
            geo = geometry.Cube8Node(pts, cubes);
            this = this@muscle.AModelConfig(geo.toCube27Node,geo);
            
            this.Version = version;
            switch this.Version    
            case 3
                types = [0 .2 .4 .6 .8 1];
                ftw = zeros(geo.GaussPointsPerElem,length(types),geo.NumElements);
                % Test: Use only slow-twitch muscles
                ftw(:,1,:) = 1;
                this.FibreTypeWeights = ftw;
                this.FibreTypes = types;
            end
        end
        
        function configureModel(this, model)
            sys = model.System;
            f = sys.f;
            model.T = 1;
            model.dt = .05;
            switch this.Version
            case 1
                f.alpha = 0;
            case 4
                model.T = 6;
                model.dt = .1;
                f.alpha = 0;
                sys.ApplyVelocityBCUntil = 2;
            end
        end
        
    end
    
    methods(Access=protected)
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            % Always fix back side
            displ_dir(:,geo.Elements(1,geo.MasterFaces(4,:))) = true;
            switch this.Version
            
            case {2,3}
                %displ_dir(:,geo.Elements(1,[6:8 11 12 18:20])) = true;
            
                % Only fix the right back corner nodes in yz-direction
    %             displ_dir([2 3],geo.Elements(1,[9 18 27])) = true;
                displ_dir(:,geo.Elements(1,[9 18 27])) = true;
            end
            
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            % Determines the dirichlet velocities.
            %
            % The unit for the applied quantities is [mm/ms] = [m/s]
            geo = this.PosFE.Geometry;
            switch this.Version
            
            case 4
                % Move the whole front
                velo_dir(2,geo.Elements(1,geo.MasterFaces(3,:))) = true;
                velo_dir_val(velo_dir) = -.3;
            end
        end
        
        function anull = seta0(this, anull)
            switch this.Version
            case {2,3}
                anull(1,:,:) = 1;
            case 4
                anull(2,:,:) = 1;
            end
        end
    end
    
    methods(Static)
        function test_DebugConfigV4
            m = muscle.Model(muscle.DebugConfig(4));
            mu = 1;
            [t,y] = m.simulate(mu);
            pdf = m.getResidualForces(t,y);
            m.plot(t,y,'PDF',pdf,'Velo',true);
        end
    end
    
end

