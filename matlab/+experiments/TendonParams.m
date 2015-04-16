classdef TendonParams < experiments.AExperimentModelConfig
    
    properties(Constant)
        TendonLength = 12; % [mm]
    end
    
    properties(SetAccess=private)
        % Assume a maximal stretch of 5%
        ExperimentalStretchPercent = linspace(0,.05,10);
        ExperimentalStretchMillimeters;
        PositioningTime = 50;
        RelaxTime = 25;
    end
    
    methods
        function this = TendonParams(varargin)
            this = this@experiments.AExperimentModelConfig(varargin{:});
            this.init;
            
            this.ExperimentalStretchMillimeters = experiments.TendonParams.TendonLength...
                *this.ExperimentalStretchPercent;
            
            this.NumOutputs = 1;
            this.NumConfigurations = length(this.ExperimentalStretchMillimeters);
            
            % SEC function from GM5 muscle.
            % Assumed a total tendon/elastic component length of 45mm to
            % attain the max isometric force of 11.2N at 4.5% overall
            % stretch (model fitted to apneurosis and tendon complex,
            % we have plain tendon material here.
            % This needs to be replaced with true experimental data)
            s = tools.SiebertTendonFun(45);
            f = s.getFunction;
            lambda = 1+this.ExperimentalStretchPercent;
            this.TargetOutputValues = f(lambda);
            
            this.VelocityBCTimeFun = tools.ConstantUntil(this.PositioningTime,.01);
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            os = m.ODESolver;
            %os.RelTol = 1e-4;
            %os.AbsTol = .001;
            
            m.T = this.PositioningTime+this.RelaxTime;
            m.dt = m.T / 300;
            
            m.DefaultMu(13) = .380; % [MPa]
            m.DefaultMu(14) = 2.05; % lambda_opt, educated guess from data
            
            % No activation here
            m.DefaultMu(2) = 0;

            % IMPORTANT! Leave off, as the cache only recognized
            % parameter/input changes but not different BCs!
            m.EnableTrajectoryCaching = false;
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
            % Runs the series of isometric tests for different parameter
            % variations to determine tendon anisotropic law components
            
            %% Test on initial values - LEAD 33000 traj
%             br = logspace(6,8,15);
%             dr = linspace(6,30,15);
%             mm = logspace(5,6,15);
%             idx = [7 8 15];
%             prefix = 'b_d_mm_detailed';
%             range = Utils.createCombinations(br,dr,mm);
            
            %% Fit for CubicToLinear law
            idx = 7;
            range = linspace(1.015,1.022,12);
            prefix = 'lambda_0_fit';
            
            %% TEST PART
            c = experiments.TendonParams('Tag',prefix);
            cap = 'Test for various b,d tendon params';      
            
            m = muscle.Model(c);
            e = tools.ExperimentRunner(m);
            nmu = size(range,2);
            
            mus = repmat(m.DefaultMu,1,nmu);
            mus(idx,:) = range;
            data = e.runExperimentsCached(mus);
            
            sp = c.ExperimentalStretchMillimeters;
            %force = squeeze(data.o(:,1,:));
            force = data.o(:,:,1)';
            
            maxbest = 5;
            %maxbest = nmu;
            nbest = min(nmu,maxbest);
            diff = Norm.L2(force-repmat(c.TargetOutputValues',1,nmu))/norm(c.TargetOutputValues);
            [diff, idx] = sort(diff);
            mus = mus(:,idx);
            for k=1:nbest
                %fprintf('%d. mu(%d): b=%10.3g, d=%10.4g, M=%10.4g with %10.4g relative L^2 error\n',k,idx(k),mus(7,k),mus(8,k),mus(15,k),diff(k));
                fprintf('%d. mu(%d): lam0=%10.5g, M=%10.5g, with %10.4g relative L^2 error\n',k,idx(k),mus(7,k),mus(8,k),diff(k));
            end
            
            pm = PlotManager;
            pm.LeaveOpen = true;
            pm.UseFileTypeFolders = false;
            ax = pm.nextPlot(['tendonparam_' c.Options.Tag],...
                ['Experiment: ' cap],'stretch [mm]','forces [mN]');
            plot(ax,sp,force(:,idx(1:nbest)));
            hold(ax,'on');
            plot(ax,sp,c.TargetOutputValues,'rx');
            legend([sprintfc('%2.2g',diff(1:nbest)) {'Data'}],'Location','NorthWest');
            
            % Set label
            % compute percent
            w = c.PosFE.Geometry.Depth;
            sp_perc = (((w+sp)./w)-1)*100;
            lab = sprintf('%2.2gmm(%2.2g%%)|',reshape([sp; sp_perc],1,[]));
            set(ax,'XTick',sp,'XTickLabel', lab(1:end-1));
            
            disp('Computation times [s]:');
            disp(data.ct);
            
            pm.savePlots(c.ImgDir,'Format',{'eps','jpg'});
        end
    end
    
end

