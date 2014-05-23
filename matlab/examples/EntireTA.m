classdef EntireTA < muscle.AModelConfig
    
    methods
        function this = EntireTA
            % Single cube with same config as reference element
            s = load(fullfile(fileparts(which(mfilename)),'..','CMISS','EntireTA.mat'));
            this = this@muscle.AModelConfig(s.geo27,s.geo8);
            geo = s.geo27;
            
            %% Muscle fibre weights
            types = [0 .2 .4 .6 .8 1];
            ftw = zeros(this.PosFE.GaussPointsPerElem,length(types),geo.NumElements);
            % Test: Use only slow-twitch muscles
            ftw(:,1,:) = 1;
            this.FibreTypeWeights = ftw;
            this.FibreTypes = types;
        end
        
        function configureModel(this, model)
            % Overload this method to set model-specific quantities like
            % simulation time etc
            
            model.T = 50;
            model.dt = 1;
            f = model.System.f;
            f.alpha = .3;
            model.System.Viscosity = .1;
            os = model.ODESolver;
            os.RelTol = .001;
            os.AbsTol = .05;
        end
    end
    
    methods(Access=protected)
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            displ_dir(:,geo.Elements(6,1:9)) = true;
            displ_dir(:,geo.Elements(8,19:27)) = true;
        end
                
        function anull = seta0(this, anull)
            % One side
            anull(1,:,:) = 1;
%             anull(1,:,[1:4 8]) = 1;
%             anull(2,:,[1:4 8]) = -1;
%             
%             % Other side
%             anull(1,:,[5 7 9:11]) = 1;
%             anull(2,:,[5 7 9:11]) = 1;
%             
            anull(1,:,[6 12]) = 0;
            anull(3,:,[6 12]) = 1;
        end
    end
    
end

