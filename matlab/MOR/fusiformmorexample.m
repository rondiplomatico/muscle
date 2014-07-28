classdef fusiformmorexample < muscle.AModelConfig
    
    properties(Access=private)
        NumParts;
    end
    
    methods
        function this = fusiformmorexample
            np = 4;
            fusiformmorexample = Belly.getBelly(np,10,1,.5,2);
            this = this@muscle.AModelConfig(fusiformmorexample);
            this.NumParts = np;
            
            %% Muscle fibre weights
            types = [0 .2 .4 .6 .8 1];
            ftw = zeros(this.PosFE.GaussPointsPerElem,length(types),this.PosFE.Geometry.NumElements);
            % Test: Use only slow-twitch muscles
            ftw(:,1,:) = .4;
            ftw(:,2,:) = .05;
            ftw(:,3,:) = .05;
            ftw(:,4,:) = .1;
            ftw(:,5,:) = .2;
            ftw(:,6,:) = .2;
            this.FibreTypeWeights = ftw;
            p = models.motorunit.Pool;
            p.FibreTypes = types;
            this.Pool = p;
        end
        
        function configureModel(~, model)
            model.T = 150;
            model.dt = .1;
            f = model.System.f;
            f.alpha = @(t)0;
        end
    end
    
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            for k = geo.NumElements-3:geo.NumElements
                displ_dir(:,geo.Elements(k,geo.MasterFaces(4,:))) = true;
            end
            for k = 1:4
                displ_dir(:,geo.Elements(k,geo.MasterFaces(3,:))) = true;
            end
        end
        
        function anull = seta0(~, anull)
            anull(2,:,:) = 1;
        end
    end
   
    
end

