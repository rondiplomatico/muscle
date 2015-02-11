classdef EntireTA < muscle.AModelConfig
    
    methods
        function this = EntireTA(varargin)
            % Single cube with same config as reference element
            this = this@muscle.AModelConfig(varargin{:});
            this.init;

            %% Muscle fibre weights
%             geo = s.geo27;
%             types = [0 .2 .4 .6 .8 1];
%             ftw = zeros(this.PosFE.GaussPointsPerElem,length(types),geo.NumElements);
%             % Test: Use only slow-twitch muscles
%             ftw(:,1,:) = .4;
%             ftw(:,2,:) = .05;
%             ftw(:,3,:) = .05;
%             ftw(:,4,:) = .1;
%             ftw(:,5,:) = .2;
%             ftw(:,6,:) = .2;
%             this.FibreTypeWeights = ftw;
%             p = models.motorunit.Pool;
%             p.FibreTypes = types;
%             this.Pool = p;
        end
        
        function configureModel(this, m)
            % Overload this method to set model-specific quantities like
            % simulation time etc
            m.T = 200;
            m.dt = 1;
            m.DefaultMu(1) = .1;
            m.DefaultMu(2) = 50;
%             m.DefaultMu(4) = 4;
            os = m.ODESolver;
            os.RelTol = .01;
            os.AbsTol = .08;
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            s = load(fullfile(fileparts(which(mfilename)),'..','CMISS','EntireTA.mat'));
            geo = s.geo27;
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            % Front face
%             displ_dir(:,geo.Elements(6,geo.MasterFaces(3,:))) = true;
            % Back faces
            displ_dir(:,geo.Elements(8,geo.MasterFaces(4,:))) = true;
            displ_dir(:,geo.Elements(8,geo.MasterFaces(2,:))) = true;
        end
                
        function anull = seta0(~, anull)
            % The elements are aligned so that the fibres go in x-direction
            anull(1,:,:) = 1;
        end
    end
    
    methods(Static)
        function test_EntireTA
            m = muscle.Model(EntireTA);
            m.simulateAndPlot;
        end
    end
    
end

