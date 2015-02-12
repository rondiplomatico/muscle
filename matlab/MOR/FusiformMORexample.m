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
        
        function this = FusiformMORexample(varargin)
            % specify as a subclass of AModelConfig
            this = this@muscle.AModelConfig(varargin{:});
            % initialise -- calls getGeometry()
            this.init;
            % specify the definition of bc    
            this.NeumannCoordinateSystem = 'local'; % same as default (not necessary)
        end
        
        function configureModel(this, model)
            % change default tolerances (necessary for convergence!!)
            model.ODESolver.RelTol = 1e-3;
            model.ODESolver.AbsTol = 1e-3;
            
            % simulation time and step size
            model.T = 350;
            model.dt = 0.1;
            
            % for offline phases
            model.EnableTrajectoryCaching = true;
            
            model.Data.useFileTrajectoryData;
            model.ComputeTrajectoryFxiData = true;
            
            % change default model parameters mu
            model.DefaultMu(1) = 1;             % viscosity
            model.DefaultMu(2) = 50;            % activation duration
            model.DefaultMu(3) = 10;            % NeumannBC (max force)
            model.DefaultMu(9) = 6.352e-10;     % Mooney-Rivlin c10
            model.DefaultMu(10) = 3.627;        % Mooney-Rivlin c01
            model.DefaultMu(13) = 73;           % muscle fibre maximal force
            
            % default input index
            model.DefaultInput = 1;             % index for u (is a cell)
            
            % specify model parameters - for trajectory computation
            sys = model.System;
            % viscosity
            sys.Params(1).Range = [1e-3 10];
            sys.Params(1).Desired = 5;
            sys.Params(1).Spacing = 'log';
            % activation-/alpha- ramp
            sys.Params(2).Range = [1 300];
            sys.Params(2).Desired = 5;
            sys.Params(2).Spacing = 'log';
            % Neumann BC
            sys.Params(3).Range = [0.1 1e3];
            sys.Params(3).Desired = 5;
            sys.Params(3).Spacing = 'log';
            
            % simulation: starting at 20 ms, fully activate in mu(2) ms
            % in mu(2) ms, up to maximal value ActivationRampMax, starting at ActivationRampOffset ms
            this.ActivationRampMax = 1; % default
            this.ActivationRampOffset = 20;            
        end
        
        % Neumann BC (mu(3)) 
        % apply a force within 10 ms, up to maxvalue mu(3), starting at 0 ms
        function u = getInputs(this)
            % Returns the inputs `u(t)` of the model.
            %
            % if neumann boundary conditions are used, this input is
            % multiplied with the mu(3) parameter, which determines the
            % maximum force that is applied. u(t) determines its temporal
            % strength.
            ramp = tools.Ramp(10,1,0);
            u{1} = ramp.getFunction;
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            % Determines the Neumann forces on the boundary (matrix B)
            %
            % The unit for the applied quantities is kiloPascal [kPa]
            %
            % See also: NeumannCoordinateSystem
            geo = this.PosFE.Geometry;
            P = [];
            if any(elemidx == [geo.NumElements-3:geo.NumElements]) && faceidx == 4
                P = -1;
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
        
        function geo = getGeometry(this)
            % geometry for model configuration
            np = 12;
            geo = Belly.getBelly(np,50,3,.5,15); 
            this.NumParts = np;            
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
