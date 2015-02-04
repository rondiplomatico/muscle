classdef IsometricActivation < experiments.AExperimentModelConfig
    % Implements the isometric contraction experiment
    
    properties(Constant)
        OutputDir = fullfile(fileparts(which(mfilename('class'))),'isometricactivation');
        ImgDir = fullfile(fileparts(which(mfilename('class'))),'isometricactivation','img');
        ExperimentalStretchMillimeters = [0 1 -1 2 -2 3 -4 4 -5 5 -6];
    end
    
    properties
        ActivationTime = 100; %[ms]
        PositioningTime = 50; %[ms]
        RelaxTime = 10; %[ms]
    end
    
    properties(SetAccess=private)
        GeoNr = 1;
    end
    
    properties(Access=private)
        % Flag for different internal geo/bc configs
        internalconfig;
    end
    
    methods
        function this = IsometricActivation(geonr, internalconfig)
            Utils.ensureDir(experiments.IsometricActivation.OutputDir);
            if nargin < 2
                internalconfig = 1;
                if nargin < 1
                    geonr = 1;
                end
            end
            this = this@experiments.AExperimentModelConfig;
            this.GeoNr = geonr;
            this.internalconfig = internalconfig;
            this.NumOutputs = 2;
            this.NumConfigurations = 11;
            % Data from GM5 muscle, originally in [N], here transferred to
            % [mN]
            this.TargetOutputValues = [9.6441	0.0217
                10.7124	0.0578
                8.8006	-0.0112
                11.1749	0.1151
                7.1639	-0.0428
                11.295	0.1533
                3.5734	-0.0604
                10.8832	0.1969
                2.1067	-0.061
                10.3544	0.2547
                1.065	-0.0568]*1000;
            
            this.VelocityBCTimeFun = tools.ConstantUntil(this.PositioningTime,.01);
        end
        
        function configureModel(this, m)
            os = m.ODESolver;
            os.RelTol = 1e-4;
            os.AbsTol = .001;
            
            m.T = this.PositioningTime + this.RelaxTime + this.ActivationTime;
            m.dt = m.T / 300;
            
            m.DefaultMu(13) = 300; % [kPa]
            m.DefaultMu(14) = 1.1; % lambda_opt, educated guess from data
            
            % Set to activation within 10ms
            m.DefaultMu(2) = 10;
            % The f function calls getAlphaRamp with mu(2), this will set
            % the default offset time correctly
            this.ActivationRampOffset = this.PositioningTime + this.RelaxTime;

            % IMPORTANT! Leave off, as the cache only recognized
            % parameter/input changes but not different BCs!
            m.EnableTrajectoryCaching = false;
        end
        
        function o = getOutputOfInterest(this, t, y)
            m = this.Model;
            % Get the residual dirichlet forces
            df = m.getResidualForces(t,y);
            
%             m.plot(t,y,'DF',df,'F',4);
            % Get the position at which to determine the passive forces
            passive_pos = floor((this.PositioningTime + this.RelaxTime)/m.dt);
            % Get the range at which to determine the max forces
            act_range = passive_pos + 1 : size(y,2);
            switch this.GeoNr
                case 1
                    idx = m.getVelocityDirichletBCFaceIdx(3,2);
                    o(1) = max(abs(sum(df(idx,act_range),1)));
                    o(2) = abs(sum(df(idx,passive_pos)));
                    
%                     idx = m.getPositionDirichletBCFaceIdx(1,1);
%                     o(3) = sum(df(idx,passive_pos));
%                     o(4) = max(sum(df(idx,act_range),1));
            end
%             figure;
%             plot(t,sum(df(idx,:)));
        end
    end
    
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.PosFE.Geometry;
            switch this.GeoNr
                case 1
                    switch this.internalconfig
                        case 1
                            % Fix front in x direction, center point in all
                            % directions
                            displ_dir(1,geo.Elements(1,geo.MasterFaces(1,:))) = true;
                            displ_dir(:,geo.Elements(1,geo.MasterFaces(1,5))) = true;
                        case {2, 3}
                            % Fix left side
                            displ_dir(:,geo.Elements(1,geo.MasterFaces(1,:))) = true;
                            % Fix right side in yz
                            displ_dir([2 3],geo.Elements(3,geo.MasterFaces(2,:))) = true;
                    end
                case 2
            end
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            geo = this.PosFE.Geometry;
            
            % Get difference in total width using the total stretch per configuration
            velo = this.ExperimentalStretchMillimeters(this.CurrentConfigNr)/this.PositioningTime;
            switch this.GeoNr
                case 1
                    velo_dir(1,geo.Elements(3,geo.MasterFaces(2,:))) = true;
                    velo_dir_val(velo_dir) = velo;
                case 2
            end
        end
        
        function anull = seta0(this, anull)
            switch this.GeoNr
                case 1
                    if this.internalconfig == 3
                        % Fibres in xz direction
                        anull([1 3],:,:) = 1;
                    else
                        % Fibres in x direction
                        anull(1,:,:) = 1;
                    end
                case 2
            end
        end
    end
    
    methods(Static)
        function runExperiments
            fi = fullfile(experiments.IsometricActivation.OutputDir,'model.mat');
            if exist(fi,'file') == 2
                load(fi);
            else
                c = experiments.IsometricActivation;
                m = muscle.Model(c);
                e = experiments.ExperimentRunner(m);
                o = e.runExperiment;
                save(fi,'m','o','c');
            end
            
            sp = experiments.IsometricActivation.ExperimentalStretchMillimeters;
            [sp,idx] = sort(sp);
            pm = PlotManager;
            pm.LeaveOpen = true;
            ax = pm.nextPlot('isomet_res','Isometric simulation results','stretch','forces [mN]');
            plot(ax,sp,o(idx,1)-o(idx,2),'r',sp,o(idx,2),'b');
            hold(ax,'on');
            tv = c.TargetOutputValues(idx,:);
            plot(ax,sp,tv(:,1)-tv(:,2),'rx',sp,tv(:,2),'bx');
            legend('Active force','Passive force','Experiment','Experiment');
            %% Set label
            % compute percent
            w = c.PosFE.Geometry.Width;
            sp_perc = (((w+sp)./w)-1)*100;
            lab = sprintf('%dmm (%g%%)|',reshape([sp; sp_perc],1,[]));
            set(ax,'XTick',sp,'XTickLabel', lab(1:end-1));
            
            pm.done;
        end
        
        function runPMAXExperiments(version)
            % Runs the series of isometric tests for different Pmax values.
            % Several variants can be chosen:
            %
            %% Geometry 1:
            % Version=1: Straight fibres in x direction, "loose" ends that
            % allow the geometry to expand when compressed
            % Version=2: Straight fibres in x direction, but completely
            % fixed ends to ensure comparability with version 3
            % Version=3: Diagonal fibres in xz direction with completely
            % fixed ends
            %
            %% Geometry 2:
            % -
            if nargin < 1
                version = 1;
            end
            
            switch version
                case 1
                    % "Normal" fibre direction only x
                    name = 'pmax';
                    cap = 'Movable setup with x-aligned fibres';
                case 2
                    % "Normal" fibre direction only x
                    name = 'pmax_fixed';
                    cap = 'Fixed setup with x-aligned fibres';
                case 3
                    % Requires the fibre direction to be xz
                    name = 'pmax_fixed_xzfibre';
                    cap = 'Fixed setup with xz-diagonal fibres';
            end
            fi = fullfile(experiments.IsometricActivation.OutputDir,['model_' name '.mat']);
            if exist(fi,'file') == 2
                load(fi);
            else
                c = experiments.IsometricActivation(1,version);
                m = muscle.Model(c);
                e = experiments.ExperimentRunner(m);
                mus = repmat(m.DefaultMu,1,4);
                mus(13,:) = 300:50:450;
                o = e.runExperiments(mus);
                save(fi,'m','o','c','mus');
            end
            
            sp = experiments.IsometricActivation.ExperimentalStretchMillimeters;
            [sp,idx] = sort(sp);
            passive = o(:,idx,2);
            active = o(:,idx,1)-passive;
            
            pm = PlotManager;
            pm.LeaveOpen = true;
            pm.UseFileTypeFolders = false;
            ax = pm.nextPlot(['isomet_' name],...
                ['Pmax' sprintf('-%d',mus(13,:)) ' experiment: ' cap],'stretch','forces [mN]');
            plot(ax,sp,active','r',sp,passive','b');
            hold(ax,'on');
            tv = c.TargetOutputValues(idx,:);
            plot(ax,sp,tv(:,1)-tv(:,2),'rx',sp,tv(:,2),'bx');
            %% Set label
            % compute percent
            w = c.PosFE.Geometry.Width;
            sp_perc = (((w+sp)./w)-1)*100;
            lab = sprintf('%dmm (%g%%)|',reshape([sp; sp_perc],1,[]));
            set(ax,'XTick',sp,'XTickLabel', lab(1:end-1));
            
            disp('Values:');
            disp(active);
            disp('Ratios:');
            disp(bsxfun(@rdivide,active,active(1,:)));
            
            pm.savePlots(experiments.IsometricActivation.ImgDir,'Format',{'eps','jpg'});
        end
    end
    
end

