classdef CubeMove < models.muscle.AMuscleConfig
    
    methods
        function this = CubeMove(varargin)
            % Single cube with same config as reference element
            this = this@models.muscle.AMuscleConfig(varargin{:});
            this.init;
            this.VelocityBCTimeFun = general.functions.Sinus(10,0,.75);
        end
        
        function configureModel(this, m)
            configureModel@models.muscle.AMuscleConfig(this, m);
            m.T = 150;
            m.dt = .5;
            m.DefaultMu(1) = .5;
            m.DefaultMu(3) = .1;
            
            % MOR pre-setup.
            % If you assign a new SpaceReducer instance, dont forget to set
            % the TargetDimensions property accordingly (is now
            % automatically set within the configureModel base routine)
            %s = spacereduction.PODGreedy;
            %s.Eps = 1e-9;
            %s.MaxSubspaceSize = 500;
            s = spacereduction.PODReducer;
            s.IncludeInitialSpace = true;
%             s.IncludeBSpan = true;
            s.Mode = 'abs';
            m.SpaceReducer = s;
        end
        
        function configureModelFinal(this)
            % Set desired reduction to full state space dimension by
            % default
            m = this.Model;
            m.SpaceReducer.Value = m.System.NumDerivativeDofs;
            configureModelFinal@models.muscle.AMuscleConfig(this);
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            if this.Options.GeoNr == 2
                geo = fem.geometry.RegularHex27Grid([0 2.5 5],-10:20,[0 2.5 5]);
            else
                geo = fem.geometry.RegularHex27Grid([0 5],-1:10,[0 5]);
            end
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.FEM.Geometry;
            % Fix all on left
            if this.Options.GeoNr == 2
                displ_dir(:,geo.Elements([59 60 119 120],geo.MasterFaces(4,:))) = true;
            else
                displ_dir(:,geo.Elements(end,geo.MasterFaces(4,:))) = true;
            end
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.FEM.Geometry;
            % Fix all on left and only the y,z directions of the back right
            % nodes
            if this.Options.GeoNr == 2
                velo_dir(2,geo.Elements([1 2 61 62],geo.MasterFaces(3,:))) = true;
            else
                velo_dir(2,geo.Elements(1,geo.MasterFaces(3,:))) = true;
            end
            velo_dir_val(velo_dir) = -.05;
        end
        
        function anull = seta0(~, anull)
           % Direction is xz
           %anull([1 3],:,:) = 1;
        end
    end
    
    methods(Static)
        function test_CubePull
            m = models.muscle.Model(models.muscle.examples.CubePull);
            mu = m.getRandomParam;
            mu(3) = .1;
            m.simulateAndPlot(true,mu,1);
        end
    end
end

