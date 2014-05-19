classdef EntireTA < muscle.AModelConfig
    
    methods
        function this = EntireTA
            % Single cube with same config as reference element
            s = load(fullfile(fileparts(which(mfilename)),'CMISS','EntireTA.mat'));
            this = this@muscle.AModelConfig(s.geo20,s.geo8);
            geo = s.geo20;
            
            %% Muscle fibre weights
            types = [0 .2 .4 .6 .8 1];
            ftw = zeros(geo.GaussPointsPerElem,length(types),geo.NumElements);
            % Test: Use only slow-twitch muscles
            ftw(:,1,:) = 1;
            this.FibreTypeWeights = ftw;
            this.FibreTypes = types;
        end
        
        function configureModel(this, model)
            % Overload this method to set model-specific quantities like
            % simulation time etc
            
            model.T = 150;
            model.dt = 1;
            f = model.System.f;
            f.alpha = .3;
            f.Viscosity = .1;
            os = model.ODESolver;
            os.RelTol = .3;
            os.AbsTol = .1;
        end
    end
    
    methods(Access=protected)
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            displ_dir(:,geo.Elements(6,1:8)) = true;
            displ_dir(:,geo.Elements(8,13:20)) = true;
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            %% Dirichlet conditions: Position (fix one side)
%             tq = this.PosFE;
            
%             for k = [1 9]
                % Quadratic
%                 velo_dir(1,tq.elems(k,[1:3 9 10 13:15])) = true;
%                 velo_dir(1,tq.elems(k,[3 10 15])) = true;
%             end
%             velo_dir_val(velo_dir) = -.1;

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

