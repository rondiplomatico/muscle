classdef MuscleTendonMixPullExperiment < experiments.AExperimentModelConfig
    % A simple experiment that pulls at one side with constant force but
    % varying muscle-tendon ratio at the fixed side.
    
    properties
        Variant;
    end
    
    methods
        function this = MuscleTendonMixPullExperiment(varargin)
            this = this@experiments.AExperimentModelConfig(varargin{:});
            this.init;
            
            this.NumConfigurations = 10;
            this.NumOutputs = 20;
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            m.T = 99;
            m.dt = 1;
            m.DefaultInput = 1;
            
            os = m.ODESolver;
            os.RelTol = 1e-5;
            os.AbsTol = 1e-5;
            
            mu = m.DefaultMu;
            mu(2) = 0;
            mu(3) = .1;
            m.DefaultMu = mu;
        end
        
        function u = getInputs(this)
            % Ramp up the external pressure
            u{1} = this.getAlphaRamp(80,1);
        end
        
        function tmr = getTendonMuscleRatio(this, points)
            % Returns the [0,1] ratio between tendon and muscle at all
            % specified points
            %
            % This method simply returns an all-zero ratio, meaning muscle only. 
            tmr = zeros(1,size(points,2));
            ratio = (this.CurrentConfigNr-1)/(this.NumConfigurations-1);
            tmr(points(1,:)<3) = ratio;
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is megaPascal [MPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            if elemidx == this.Options.GeoNr && faceidx == 2
                P = 1;
            end
        end
        
        function o = getOutputOfInterest(this, ~, y)
            m = this.Model;
            facedof = m.getFaceDofsGlobal(this.Options.GeoNr,2,1);
            % Get x position of one node at end as output
            o = NaN(1,this.NumOutputs);
            idx = floor(linspace(1,m.T,this.NumOutputs));
            hasdata = idx <= size(y,2);
            o(hasdata) = y(facedof(1),idx(hasdata));
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            switch this.Options.GeoNr
                case 1
                    % Four cubes in a row
                    [pts, cubes] = geometry.Cube8Node.DemoGrid([0 10],[0 10],[0 10]);
                case 2
                    % Four cubes in a row
                    [pts, cubes] = geometry.Cube8Node.DemoGrid([0 3 10],[0 10],[0 10]);
            end
            geo = geometry.Cube8Node(pts, cubes);
            geo = geo.toCube27Node;    
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            %% Dirichlet conditions: Position (fix one side)
            geo = this.PosFE.Geometry;
            displ_dir(1,geo.Elements(1,geo.MasterFaces(1,:))) = true;
        end
        
        function anull = seta0(~, anull)
            anull(1,:,:) = 1;
        end
    end
    
    methods(Static)
        function test_MuscleTendonMixPullExperiment(geo)
            if nargin < 1
                geo = 1;
            end
            %% Init
            c = experiments.MuscleTendonMixPullExperiment('GeoNr',geo);
            m = c.createModel;
            e = tools.ExperimentRunner(m);
            rules = [3 5 7 10 15 20];
            nrules = length(rules);
            allo = zeros(c.NumConfigurations,c.NumOutputs,nrules);
            ctimes = zeros(c.NumConfigurations,nrules);
            
            %% Computation
            for gr = 1:nrules
                grule = rules(gr);
                m.setGaussIntegrationRule(grule);
                [o, ctimes(:,gr)] = e.runExperiment;
                allo(:,:,gr) = o;
            end
            %% Save
            save(fullfile(c.OutputDir,['output_' c.getOptionStr]),'allo','ctimes','rules');
            %% Load
%             load(fullfile(c.OutputDir,['output_' c.getOptionStr '.mat']));
            
            %% Setup
            m.plotGeometrySetup;
            %% per-gauss-rule comparison
            pm = PlotManager(false,1,2);
            pm.LeaveOpen = true;
            tmr = ((1:c.NumConfigurations)-1)/(c.NumConfigurations-1);
            [TMR,T] = meshgrid(tmr,...
                    floor(linspace(1,m.T,c.NumOutputs)));  
            for gr = 1:nrules
                grule = rules(gr);
                ax = pm.nextPlot(sprintf('gaussrule_%d',grule),...
                    sprintf('Right face position over time and tm-ratio, Gauss %d-point rule',grule),...
                    'tendon-muscle ratio [0,1]','time [ms]');
                surfc(TMR,T,allo(:,:,gr)','Parent',ax,'FaceColor','interp');
                ax = pm.nextPlot(sprintf('gaussrule_%d_err',grule),...
                    sprintf('Gauss %d-point rule error compared to %d-point rule',grule,rules(end)),...
                    'tendon-muscle ratio [0,1]','time [ms]');
%                 err = allo(:,:,gr)-allo(:,:,end)
                err = abs((allo(:,:,gr)-allo(:,:,end)))./allo(:,:,end);
                mesh(TMR,T,err','Parent',ax);
            end
            %% End position plot
            pm = PlotManager;
            pm.LeaveOpen = true;
            ax = pm.nextPlot('endpos','Comparison of end-positions over TMR and gauss rules',...
                'tendon-muscle ratio [0,1]','gauss integration rule [points]');
            [TMR,GR] = meshgrid(tmr,1:nrules);
            surfc(TMR,GR,squeeze(allo(:,end,:))','Parent',ax,'FaceColor','interp');
            pm.done;
            set(ax,'YTickLabel',sprintfc('%d',rules),'YTick',1:nrules);
        end
    end
    
end

