classdef CubePull < muscle.AModelConfig
    
    methods
        function this = CubePull(varargin)
            % Single cube with same config as reference element
            this = this@muscle.AModelConfig(varargin{:});
            this.init;
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            m.T = 7;
            m.dt = .005;
            m.DefaultMu(1) = .1;
            m.DefaultMu(3) = 100;
            m.DefaultInput = 1;
            %os = m.ODESolver;
            %os.RelTol = .0001;
            %os.AbsTol = .05;
        end
        
        function P = getBoundaryPressure(~, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is kiloPascal [kPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            if elemidx == 1 && faceidx == 3
                P = -1;
            end
        end
        
        function u = getInputs(this)
            u = {this.getAlphaRamp(1,1)};
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(~)
            %[pts, cubes] = geometry.Cube8Node.DemoGrid([0 2.5 5],-1:20,[0 2.5 5]);
            [pts, cubes] = geometry.Cube8Node.DemoGrid([0 5],-1:10,[0 5]);
            geo = geometry.Cube8Node(pts, cubes);
            geo = geo.toCube27Node;
            %geo = geo.toCube20Node;
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            % Fix all on left and only the y,z directions of the back right
            % nodes
            displ_dir(2,geo.Elements(end,geo.MasterFaces(4,:))) = true;
            displ_dir(:,geo.Elements(end,geo.MasterFaces(4,1))) = true;
        end
        
        function anull = seta0(~, anull)
           % Direction is xz
           %anull([1 3],:,:) = 1;
        end
    end
    
    methods(Static)
        function m = test_CubePull_MOR
            m = muscle.Model(CubePull);
            %m.EnableTrajectoryCaching = true;
            
            %m.ComputeParallel = true;
            %m.Data.useFileTrajectoryData;
            
            forces = linspace(-200,100,5);
            %mus = repmat(m.DefaultMu,1,length(forces));
            %mus(3,:) = forces;
            
            m.Sampler = sampling.ManualSampler(forces);
            m.TrainingParams = 3;
            m.TrainingInputs = 1;
            
            s = spacereduction.PODReducer;
            s.Mode = 'abs';
            s.Value = m.System.StateSpaceDimension;
            m.SpaceReducer = s;
            
            m.offlineGenerations;
            
            A = m.Data.TrajectoryData.toMemoryMatrix;
            U = A(1:610,:);
            V = A(611:1220,:);
            [uu,us,uv] = svd(U,'econ');
            [vu,vs,vv] = svd(V,'econ');
            pm = PlotManager;
            ma.plotReductionOverview(pm);
            ax = pm.nextPlot('singvals_split',...
                'Singular values of same training data when split into u/v parts','singular value number','value');
            semilogy(ax,1:610,diag(us),'r',1:610,diag(vs),'b')
            pm.done;
            pm.savePlots(pwd,'Format',{'jpg','pdf'});
        end
    end
    
end

