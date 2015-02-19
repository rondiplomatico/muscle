classdef TendonParams < experiments.AExperimentModelConfig
    
    properties(Constant)
        ExperimentalStretchMillimeters = linspace(0,.5,7);
    end
    
    properties
        PositioningTime = 50;
        RelaxTime = 25;
    end
    
    methods
        function this = TendonParams(varargin)
            this = this@experiments.AExperimentModelConfig(varargin{:});
            this.init;
            
            this.NumOutputs = 1;
            this.NumConfigurations = length(this.ExperimentalStretchMillimeters);
            % Data from GM5 muscle, originally in [N], here transferred to
            % [mN]
            s = tools.SiebertTendonFun;
            f = s.getFunction;
            this.TargetOutputValues = f(this.ExperimentalStretchMillimeters);
            
            this.VelocityBCTimeFun = tools.ConstantUntil(this.PositioningTime,.01);
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            os = m.ODESolver;
            %os.RelTol = 1e-4;
            %os.AbsTol = .001;
            
            m.T = this.PositioningTime+this.RelaxTime;
            m.dt = m.T / 300;
            
            m.DefaultMu(13) = 380; % [kPa]
            m.DefaultMu(14) = 2.05; % lambda_opt, educated guess from data
            
            % No activation here
            m.DefaultMu(2) = 0;

            % IMPORTANT! Leave off, as the cache only recognized
            % parameter/input changes but not different BCs!
            m.EnableTrajectoryCaching = false;
            
            m.System.f.MarkertMaxModulus = 1e6; % [mN]
        end
        
        function configureModelFinal(this)
            m = this.Model;
            m.Plotter.DefaultArgs = {'Fibres',false};
        end
        
        function tmr = getTendonMuscleRatio(this, points)
            % Have tendon only here!
            tmr = ones(1,size(points,2));
        end
        
        function o = getOutputOfInterest(this, t, y)
            m = this.Model;
            % Get the residual dirichlet forces at end time
            df = m.getResidualForces(t(end),y(:,end));
            idx = m.getPositionDirichletBCFaceIdx(1:4,3);
            o(1) = sum(df(idx));
            %o(2) = sum(df(idx,passive_pos));            
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            cent = -12;
            k = kernels.GaussKernel(10);
            k2 = kernels.GaussKernel(12);
            rad = @(t)[k.evaluate(t,cent) k2.evaluate(t,cent)]';
            geo = Belly.getBelly(3,12,'Radius',rad,'InnerRadius',0.74);
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.PosFE.Geometry;
            displ_dir(2,geo.Elements(1:4,geo.MasterFaces(3,:))) = true;
            displ_dir(:,1) = true;
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            geo = this.PosFE.Geometry;
            % Get difference in total width using the total stretch per configuration
            velo = this.ExperimentalStretchMillimeters(this.CurrentConfigNr)/this.PositioningTime;
            velo_dir(2,geo.Elements(9:12,geo.MasterFaces(4,:))) = true;
            velo_dir_val(velo_dir) = velo;
            
            % Add slight movement in x direction for numerical treatment of
            % compression cases
%             extra = false(size(velo_dir));
%             extra(1,geo.Elements(9:12,geo.MasterFaces(4,:))) = true;
%             velo_dir = velo_dir || extra;
%             velo_dir_val(extra) = 0.001;
        end
        
        function anull = seta0(this, anull)
            anull(2,:,:) = 1;
        end
    end
    
    methods(Static)
        
        function runExperiments()
            % Runs the series of isometric tests for different Pmax values.
            % Several variants can be chosen:
            %
            
            %% Initial rough test
%             br = logspace(0,6,10);
%             dr = linspace(3,40,10);
%             prefix = 'b_d_100';
            
            %% Initial rough test
            br = logspace(log10(4.64e+04),6,10);
            dr = linspace(7,30,10);
            prefix = 'b_d_refined_100';

            c = experiments.TendonParams('Tag',prefix);
            cap = 'Test for various b,d tendon params';
            
            %% NOT CONFIGURABLE PART
            range = Utils.createCombinations(br,dr);
            mustr = [sprintf('-%g',br) '/' sprintf('-%g',dr)];
            idx = [7 8];            
            
            m = muscle.Model(c);
            % Use mooney-rivlin for muscle material here
            m.DefaultMu(11:12) = m.DefaultMu(9:10);
            e = tools.ExperimentRunner(m);
%             e.RunParallel = false;
            mus = repmat(m.DefaultMu,1,size(range,2));
            mus(idx,:) = range;
            data = e.runExperimentsCached(mus);
            
            sp = c.ExperimentalStretchMillimeters;
            force = data.o(:,:,1)';
            
            nbest = 6;
            diff = Norm.L2(force-repmat(c.TargetOutputValues',1,100))/norm(c.TargetOutputValues);
            [diff, idx] = sort(diff);
            mus = mus(:,idx);
            for k=1:nbest
                fprintf('%d. b=%10.3g, d=%10.4g with %10.4g relative L^2 error\n',k,mus(7,k),mus(8,k),diff(k));
            end
            
            pm = PlotManager;
            pm.LeaveOpen = true;
            pm.UseFileTypeFolders = false;
            ax = pm.nextPlot(['tendonparam_' c.Options.Tag],...
                ['B/D' mustr ' experiment: ' cap],'stretch [mm]','forces [mN]');
            plot(ax,sp,force(:,idx(1:nbest)));
            hold(ax,'on');
            plot(ax,sp,c.TargetOutputValues,'rx');
            legend([sprintfc('%d.',1:nbest) {'Data'}]);
            
            % Set label
            % compute percent
            w = c.PosFE.Geometry.Width;
            sp_perc = (((w+sp)./w)-1)*100;
            lab = sprintf('%dmm (%g%%)|',reshape([sp; sp_perc],1,[]));
            set(ax,'XTick',sp,'XTickLabel', lab(1:end-1));
            
%             disp('Computation times [s]:');
%             disp(data.ct);
            
            pm.savePlots(c.ImgDir,'Format',{'eps','jpg'});
        end
    end
    
end

