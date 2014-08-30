classdef LongMORexample < muscle.AModelConfig
% A long geometry with 20% deviation from default cubic positions and
% complex fibre structure

    properties(Constant)
       OutputDir = fullfile(fileparts(which(mfilename)),'LongMOR_output') 
    end
    %
    %%
    methods
        function this = LongMORexample(devi)
            if nargin < 1
                devi = .2;
            end
            % Single cube with same config as reference element
            [pts, cubes] = geometry.Cube20Node.DemoGrid(-10:2.5:10,-40:2.5:40, 0:2.5:10, devi);
            geo = geometry.Cube20Node(pts, cubes);
            this = this@muscle.AModelConfig(geo);
        end
        
        function configureModel(this, model)
            %
            model.T = 100;
            model.dt = 1;
            model.EnableTrajectoryCaching = true;
            
            model.Data.useFileTrajectoryData;
            model.ComputeTrajectoryFxiData = true;
            
            % default model parameters 
            model.DefaultMu = [1; 50];      % mu = [viscosity; activation duration]
            
            % specify model parameters (mu = [viscosity; activation duration])
            sys = model.System;
            %viscosity
            sys.Params(1).Range = [0.5 10];
            sys.Params(1).Desired = 10;
            sys.Params(1).Spacing = 'log';
            % activation duration
            sys.Params(2).Name = 'alpha-ramp';
            sys.Params(2).Range = [20 100];
            sys.Params(2).Desired = 10;
            sys.Params(2).Spacing = 'log';
            
            f = sys.f;
            % Material set (see main comment)
            f.c10 = 6.352e-10; % [kPa]
            f.c01 = 3.627; % [kPa]
            f.b1 = 2.756e-5; % [kPa]
            f.d1 = 43.373; % [-]

            f.Pmax = 73; % [kPa]
            f.lambdafopt = 1.2; % [-]            
        end
        
        function prepareSimulation(this, mu, inputidx)
            % 
            sys = this.Model.System;
            sys.f.alpha = this.getAlphaRamp(mu(2),1);    % (in ..ms, up to maxvalue.., starting at ..ms)
        end
        
        function o = getOutputOfInterest(this, model, t, uvw)
            % Writes the data of interest into o
            %
            geo = model.Config.PosFE.Geometry;
            %[df,nf] = model.getResidualForces(t,uvw);
            uvw = model.System.includeDirichletValues(t,uvw);
            % get displacement of loose/right end
            facenode_idx = [];
            for k = geo.NumElements-3:geo.NumElements
                facenode_idx = [facenode_idx; model.getFaceDofsGlobal(k,4,2)];
            end
            facenode_idx = unique(facenode_idx);
%             o = uvw(facenode_idx,:);            % gives matrix, where each rowvector shows positon of one node over time
            o = mean(uvw(facenode_idx,:),1);      % gives rowvector of mean node-position for all timesteps
        end
        
    end
    %
    %%
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            % Dirichlet conditions: Position (fix both sides)
            geo = this.PosFE.Geometry;
            %left
            for k = 1:32
                displ_dir(:,geo.Elements(k,geo.MasterFaces(3,:))) = true;
            end
            %right
            for k = geo.NumElements-31:geo.NumElements
                displ_dir(:,geo.Elements(k,geo.MasterFaces(4,:))) = true;
            end
        end
        
        function anull = seta0(this, anull)
            fe = this.PosFE;
            if fe.GaussPointsPerElem ~= 27
                warning('a0 designed for 27 gauss points!');
            end
            x = linspace(0,1,8*3);
            basea0 = [sin(x*pi); cos(x*pi)]; %; zeros(size(x))
            front = fe.GaussPoints(2,:) < 0;
            mid = fe.GaussPoints(2,:) == 0;
            back = fe.GaussPoints(2,:) > 0;
            % Direction is x
            for m = 1:8
                off = (m-1)*3;
                anull([1 3],front,[m m+8]) = basea0(2,off+1);
                anull([1 3],mid,[m m+8]) = basea0(2,off+2);
                anull([1 3],back,[m m+8]) = basea0(2,off+3);
                
                anull(2,front,[m m+8]) = basea0(1,off+1);
                anull(2,mid,[m m+8]) = basea0(1,off+2);
                anull(2,back,[m m+8]) = basea0(1,off+3);
            end
        end
        
    end
    
end