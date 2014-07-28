classdef FusiformMORexample < muscle.AModelConfig
    
    properties(Access=private)
        % number of parts in one quarter of the belly
        NumParts;
        
        Loads;
        Pressures;
    end
    
    methods
        function this = FusiformMORexample
            % setting up the fusiform muscle geometry
            np = 4;
            belly = Belly.getBelly(np,10,1,.5,2);
            this = this@muscle.AModelConfig(belly);
            this.NumParts = np;
            % default model parameters
            this.Model.DefaultMu = [1; 1];
            
            % specify the definition of bc    
            this.NeumannCoordinateSystem = 'local';
        end
        
        function prepareSimulation(this, mu, inputidx)
            % 
            f = this.Model.System.f;
            f.alpha = this.getAlphaRamp(mu(2),1);
            
        end
        
        function configureModel(this, model)
            model.T = 150;
            model.dt = 0.1;
            
            f = model.System.f;
            % Material set (see main comment)
            f.c10 = 6.352e-10; % [kPa]
            f.c01 = 3.627; % [kPa]
            f.b1 = 2.756e-5; % [kPa]
            f.d1 = 43.373; % [-]

            f.Pmax = 250; % [kPa]
            f.lambdafopt = 1; % [-]
           
            os = model.ODESolver;
            os.RelTol = 0.01;
            os.AbsTol = 0.1;
            
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            % Determines the Neumann forces on the boundary (matrix B)
            %
            % The unit for the applied quantities is kiloPascal [kPa]
            %
            % See also: NeumannCoordinateSystem
            P = [];
            if any(elemidx == [13:16]) && faceidx == 4
                P = 1;
            end
        end
        
        function u = getInputs(this)
            % vector u in nl dynamical system
            u = {@(t)1};
        end
        
    end
    
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
   
    % Test
    methods(Static)
        
        function runTest
            modconf = FusiformMORexample;
            geo = modconf.PosFE.Geometry;
            model = muscle.Model(modconf);
            [t,y] = model.simulate;
        end

    
    end

    
    
end

