classdef IsometricActivation < experiments.AExperimentModelConfig
    % Implements the isometric contraction experiment
    
    properties(Constant)
        OutputDir = fullfile(fileparts(which(mfilename('class'))),'isometricactivation');
        ExperimentalStretchPercents = [0 1 -1 2 -2 3 -4 4 -5 5 -6];
    end
    
    properties
        ActivationTime = 100; %[ms]
        PositioningTime = 50; %[ms]
        RelaxTime = 10; %[ms]
    end
    
    properties(SetAccess=private)
        GeoNr = 1;
    end
    
    methods
        function this = IsometricActivation
            Utils.ensureDir(experiments.IsometricActivation.OutputDir);
            this = this@experiments.AExperimentModelConfig;
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
        end
        
        function configureModel(this, m)
            os = m.ODESolver;
            os.RelTol = 1e-4;
            os.AbsTol = .001;
            
            m.T = this.PositioningTime + this.RelaxTime + this.ActivationTime;
            m.dt = m.T / 300;
            
            f = m.System.f;
            f.Pmax = 1250; % [kPa], in Paper 25N/cm², but kPa = 0.1N/cm² 
            f.lambdafopt = 1.3; % [-]
            % Set to activation within 10ms
            m.DefaultMu(2) = 10;
            % The f function calls getAlphaRamp with mu(2), this will set
            % the default offset time correctly
            this.ActivationRampOffset = this.PositioningTime + this.RelaxTime;
            
            m.System.ApplyVelocityBCUntil = this.PositioningTime;

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
                    idx = m.getVelocityDirichletBCFaceIdx(1,2);
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
                    % Fix front
                    displ_dir(:,geo.Elements(1,geo.MasterFaces(1,:))) = true;
                case 2
            end
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            geo = this.PosFE.Geometry;
            
            % Get difference in total width using the total stretch per configuration
            deltawidth = geo.Width*this.ExperimentalStretchPercents(this.CurrentConfigNr)/100;
            switch this.GeoNr
                case 1
                    velo_dir(1,geo.Elements(1,geo.MasterFaces(2,:))) = true;
                    velo_dir_val(velo_dir) = deltawidth/this.PositioningTime;
                case 2
            end
        end
        
        function anull = seta0(~, anull)
            % Fibres in x direction
            anull(1,:,:) = 1;
        end
    end
    
    methods(Static)
        function runExperiments
            fi = fullfile(experiments.IsometricActivation.OutputDir,'model.mat');
            loaded = false;
            if exist(fi,'file') == 2
                load(fi);
                loaded = true;
                c = m.Config;
            else
                c = experiments.IsometricActivation;
                m = muscle.Model(c);
                e = experiments.ExperimentRunner(m);
                o = e.runExperiment;
            end
            
            sp = experiments.IsometricActivation.ExperimentalStretchPercents;
            [sp,idx] = sort(sp);
            pm = PlotManager;
            pm.LeaveOpen = true;
            ax = pm.nextPlot('isomet_res','Isometric simulation results','stretch percent [%]','forces [mN]');
            plot(ax,sp,o(idx,1)-o(idx,2),'r',sp,o(idx,2),'b');
            hold(ax,'on');
            tv = c.TargetOutputValues(idx,:);
            plot(ax,sp,tv(:,1)-tv(:,2),'rx',sp,tv(:,2),'bx');
            legend('Active force','Passive force','Experiment','Experiment');
            pm.done;
            
%             if ~loaded
%                 save(fi,'m','o');
%             end
        end
    end
    
end

