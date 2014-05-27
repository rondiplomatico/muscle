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
            ftw(:,1,:) = .4;
            ftw(:,2,:) = .05;
            ftw(:,3,:) = .05;
            ftw(:,4,:) = .1;
            ftw(:,5,:) = .2;
            ftw(:,6,:) = .2;
            this.FibreTypeWeights = ftw;
            p = motorunit.Pool;
            p.FibreTypes = types;
            this.Pool = p;
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
%             displ_dir(:,geo.Elements(6,geo.MasterFaces(3,:))) = true;
            displ_dir(:,geo.Elements(8,geo.MasterFaces(4,:))) = true;
            displ_dir(:,geo.Elements(8,geo.MasterFaces(2,:))) = true;
        end
                
        function anull = seta0(this, anull)
            % One side
%             anull(1,:,:) = 1;
%             anull(2,:,:) = -1;

            anull(1,:,[1:4 8]) = 2;
            anull(2,:,[1:4 8]) = 1;
% %             
% %             % Other side
            anull(1,:,[5 7 9:11]) = 1;
            anull(2,:,[5 7 9:11]) = -2;
% %             
%             anull(1,:,[6 12]) = 0;
            anull(2,:,[6 12]) = 1;
        end
    end
    
end

