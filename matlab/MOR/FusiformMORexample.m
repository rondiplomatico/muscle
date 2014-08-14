classdef FusiformMORexample < muscle.AModelConfig
    
    properties(Constant)
       OutputDir = fullfile(fileparts(which(mfilename)),'FusiformMOR_output') 
    end
    %
    %%
    properties(Access=private)
        
        % number of parts in one quarter of the belly
        NumParts;
        
        Loads;
        Pressures;
    end
    %
    %%
    methods
        
        function this = FusiformMORexample
            % setting up the fusiform muscle geometry
%             np = 4;
%             belly = Belly.getBelly(np,10,1,.5,2);
            np = 12;
            belly = Belly.getBelly(np,50,3,.5,15);
            this = this@muscle.AModelConfig(belly);
            this.NumParts = np;
            % specify the definition of bc    
            this.NeumannCoordinateSystem = 'local';
        end
        
        function configureModel(this, model)
            model.T = 250;
            model.dt = 0.1;
            model.EnableTrajectoryCaching = true;
            
            model.Data.useFileTrajectoryData;
            model.ComputeTrajectoryFxiData = true;
            
            % default model parameters 
            model.DefaultMu = [1; 50; -10];      % mu = [viscosity; activation duration; NeumannBC (max force)]
            model.DefaultInput = 1;             % index for u (is a cell)
            
            % specify model parameters (mu = [viscosity; activation duration; NeumannBC (max force)])
            sys = model.System;
            sys.Params(1).Range = [1e-3 5];
            sys.Params(1).Desired = 5;
            sys.Params(2).Name = 'alpha-ramp';
            sys.Params(2).Range = [1 200];
            sys.Params(2).Desired = 5;
            sys.addParam('Neumann BC', [-1e3 0], 5);
            
            f = model.System.f;
            % Material set (see main comment)
            f.c10 = 6.352e-10; % [kPa]
            f.c01 = 3.627; % [kPa]
            f.b1 = 2.756e-5; % [kPa]
            f.d1 = 43.373; % [-]

            f.Pmax = 73; % [kPa]
            f.lambdafopt = 1.2; % [-]
           
            %os = model.ODESolver;
            %os.RelTol = 0.02;%0.01;
            %os.AbsTol = 0.1;%0.1;

        end
        
        function prepareSimulation(this, mu, inputidx)
            % 
            sys = this.Model.System;
            sys.f.alpha = this.getAlphaRamp(mu(2),1,20);    % (in ..ms, up to maxvalue.., starting at ..ms)
            %sys.f.alpha = @(t)0;
            sys.Inputs{1} = this.getAlphaRamp(10,mu(3));    % (in ..ms, up to maxvalue.., starting at ..ms)
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            % Determines the Neumann forces on the boundary (matrix B)
            %
            % The unit for the applied quantities is kiloPascal [kPa]
            %
            % See also: NeumannCoordinateSystem
            geo = this.PosFE.Geometry;
            P = [];
%             if any(elemidx == [13:16]) && faceidx == 4
            if any(elemidx == [geo.NumElements-3:geo.NumElements]) && faceidx == 4
                P = 1;
            end
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
            o = mean(uvw(facenode_idx,:),1);    % gives rowvector of mean node-position for all timesteps
        end
        
    end
    %
    %%
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            % Dirichlet conditions: Position (fix left side)
            geo = this.PosFE.Geometry;
            for k = 1:4
                displ_dir(:,geo.Elements(k,geo.MasterFaces(3,:))) = true;
            end
        end
        
        function anull = seta0(~, anull)
            % fibres in y direction
            anull(2,:,:) = 1;
        end
        
    end
    %
    %%
    % Testscript
    methods(Static)
        
        function runTest
            modconf = FusiformMORexample;
            geo = modconf.PosFE.Geometry;
            model = muscle.Model(modconf);
            [t,y] = model.simulate;
            model.plot(t,y);
        end

    end
    
end

