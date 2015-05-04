classdef IsometricActivation < experiments.AExperimentModelConfig
    % Implements the isometric contraction experiment
    
    properties(Constant)
        ExperimentalStretchMillimeters = [0 1 -1 2 -2 3 -4 4 -5 5 -6];
    end
    
    properties
        ActivationTime = 50; %[ms]
        PositioningTime = 50; %[ms]
        RelaxTime = 5; %[ms]
    end
    
    properties(Access=private)
        rat3;
    end
    
    methods
        function this = IsometricActivation(varargin)
            this = this@experiments.AExperimentModelConfig(varargin{:});
            this.addOption('BC',1);
            this.addOption('Activate',1);
            this.rat3 = experiments.Rat3Geo;
            this.init;
            
            this.NumOutputs = 2;
            this.NumConfigurations = length(this.ExperimentalStretchMillimeters);
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
            if this.Options.GeoNr == 4
                this.rat3.load;
            end
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            os = m.ODESolver;
            os.RelTol = 1e-4;
            os.AbsTol = .001;
            
            m.T = this.PositioningTime + this.RelaxTime;
            % Only simulate activation if enabled
            if this.Options.Activate
                m.T = m.T + this.ActivationTime;
            end
            m.dt = m.T / 300;
            
            m.DefaultMu(13) = .380; % [MPa]
            m.DefaultMu(14) = 2.05; % lambda_opt, educated guess from data
            
            % Set to activation within 10ms
            m.DefaultMu(2) = 10;
            % The f function calls getAlphaRamp with mu(2), this will set
            % the default offset time correctly
            this.ActivationRampOffset = this.PositioningTime + this.RelaxTime;

            % IMPORTANT! Leave off, as the cache only recognized
            % parameter/input changes but not different BCs!
            m.EnableTrajectoryCaching = false;
        end
        
        function configureModelFinal(this)
            if this.Options.GeoNr == 2 || this.Options.GeoNr == 3
                m = this.Model;
                m.Plotter.DefaultArgs = {'Fibres',false};
            end
        end
        
        function tmr = getTendonMuscleRatio(this, points)
            % Returns the [0,1] ratio between tendon and muscle at all
            % specified points
            %
            % This method simply returns an all-zero ratio, meaning muscle only. 
            tmr = zeros(1,size(points,2));
            o = this.Options;
            x = points(1,:);
            y = points(2,:);
            z = points(3,:);
            switch this.Options.GeoNr
                case 2
                    tmr = zeros(size(y));
                    fr = @(x)1.2*sqrt(max(0,27-x))-2;
                    right = z > fr(y);
                    tmr(right) = 1;
                    fl = @(x)-1*sqrt(x)+1;
                    left = z < fl(y);
                    tmr(left) = 1;
                case 3
                    inner = 4.5;
                    cent = [0 11 0];
                    fac = [.5 1 1];
                    fr = @(x,y,z)sqrt(fac(1)*(x-cent(1)).^2 ...
                        + fac(2)*(y-cent(2)+3*z).^2 ...
                        + fac(3)*(z-cent(3)).^2)/sqrt(3)-inner;
                    tmr = fr(x,y,z);
                    %tmr = tmr/norm(cent);
                    rate = 1;
                    tmr = min(rate,max(0,tmr))/rate;
%                     tmr
%                     right = z > fr(y);
%                     tmr(right) = 1;
%                     fl = @(x)-1*sqrt(x)+1;
%                     left = z < fl(y);
%                     tmr(left) = 1;
                case 4
                    tmr = this.rat3.getTMR(x,y,z);
            end            
        end
        
        function o = getOutputOfInterest(this, t, y)
            m = this.Model;
            
            % Get the position at which to determine the passive forces
            passive_pos = floor((this.PositioningTime + this.RelaxTime)/m.dt);
            
            % Get the residual dirichlet forces on interesting part only
            % (=passive_pos until end)
            all = passive_pos : size(y,2);
            df = m.getResidualForces(t(all),y(:,all));
            
            switch this.Options.GeoNr
                case 1
                    idx = m.getVelocityDirichletBCFaceIdx(3,2);        
                case {2,3}    
                    %this.Geometry.Nodes(2,:) > 33 & m.System.bool_v_bc_nodes
                    idx = m.getVelocityDirichletBCFaceIdx(37:48,4);
            end
            % Passive force
            o(2) = sum(df(idx,1));
            
            % Get the range at which to determine the max forces, relative
            % to passive_pos (active part begins after passiv_pos till end)
            force = sum(df(idx,2:end),1);
            % Get the force with the largest magnitude
            [~, pos] = max(abs(force));
            o(1) = force(pos);
        end
        
        function plotGeometryInfo(this, varargin)
            plotGeometryInfo@muscle.AModelConfig(this, varargin{:});
            ax = gca;
            view(ax, [90 0]);
            x = 0:.5:35;
            fr = @(x)1.2*sqrt(max(0,27-x))-2;
            plot3(ax,zeros(size(x)),x,fr(x),'LineWidth',3);
            hold(ax,'on');
            fl = @(x)-1*sqrt(x)+1;
            plot3(ax,zeros(size(x)),x,fl(x),'LineWidth',3);
            axis(ax,[-4 4 -1 36 -3 3]);
        end
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            switch this.Options.GeoNr
                case 1
                    [pts, cubes] = geometry.Cube27Node.DemoGrid(linspace(0,25,4),[0 7],[0 4]);
                    geo = geometry.Cube27Node(pts,cubes);
                case 2
                    cent = 12;
                    k = kernels.GaussKernel(10);
                    k2 = kernels.GaussKernel(12);
                    rad = @(t)[k.evaluate(t,cent)*3.5 k2.evaluate(t,cent)*2]';
                    geo = Belly.getBelly(4,35,'Radius',rad,'Layers',[.8 1],'InnerRadius',.3);
                case 3
                    cent = 12;
                    k = kernels.GaussKernel(8.5);
                    k2 = kernels.GaussKernel(9.5);
                    rad = @(t)[k.evaluate(t,cent)*3.5 k2.evaluate(t,cent)*2]';
                    geo = Belly.getBelly(4,35,'Radius',rad,'Layers',[.8 1],'InnerRadius',.3);
                case 4
                    geo = this.rat3.getGeometry;
            end
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            geo = this.PosFE.Geometry;
            o = this.Options;
            switch o.GeoNr
                case 1
                    switch o.BC
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
                case {2,3}
                    displ_dir(2,geo.Nodes(2,:) == 0) = true;
                    displ_dir(:,1) = true;
                    displ_dir([1 3],geo.Nodes(2,:) > 33) = true;
                case 4
                    displ_dir(:,geo.Elements(1:12,geo.MasterFaces(3,:))) = true;
            end
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            geo = this.PosFE.Geometry;
            
            % Get difference in total width using the total stretch per configuration
            velo = this.ExperimentalStretchMillimeters(this.CurrentConfigNr)/this.PositioningTime;
            o = this.Options;
            switch o.GeoNr
                case 1
                    velo_dir(1,geo.Elements(3,geo.MasterFaces(2,:))) = true;
                case {2,3}
                    velo_dir(2,geo.Nodes(2,:) > 33) = true;
                case 4
                    velo_dir(:,geo.Elements(49:60,geo.MasterFaces(4,:))) = true;
            end
            velo_dir_val(velo_dir) = velo;
        end
        
        function anull = seta0(this, anull)
            o = this.Options;
            switch o.GeoNr
                case 1
                    if o.BC == 3
                        % Fibres in xz direction
                        anull(1,:,:) = 1;
                        anull(3,:,:) = .3;
                    else
                        % Fibres in x direction
                        anull(1,:,:) = 1;
                    end
                case {2,3}
                    anull(2,:,:) = 1;
                    anull(3,:,:) = 1;
                case 4
                    cached = fullfile(this.OutputDir,'fibredirections_gp.mat');
                    if exist(cached,'file') == 2
                        s = load(cached);
                        anull = s.anull;
                    else
                        fe = this.PosFE;
                        g = fe.Geometry;
                        for m = 1:g.NumElements
                            gp = g.Nodes(:,g.Elements(m,:)) * fe.N(fe.GaussPoints);
                            anull(:,:,m) = this.rat3.getA0(gp);
                        end
                        save(cached,'anull');
                    end
            end
        end
    end
    
    methods(Static)
        
        function runExperiments()
            % Runs the series of isometric tests for different Pmax values.
            % Several variants can be chosen:
            %
            Geo = 4;
            
            %% pmax/lamdaopt rough test
%             pmaxr = 250:50:400;
%             lamr = 1.8:.05:2.1;
%             range = Utils.createCombinations(pmaxr,lamr);
%             idx = [13 14];
%             prefix = 'pmax_lamopt';
            
            %% Test with no anisotropic muscle force
%             range = [0; 400; 1.8];
%             idx = [5 13 14];
%             prefix = 'zero_muscle_aniso';
            
            %% Test with no anisotropic muscle force and various mr-coeffs
            pmaxr = 400;
            lamr = 2.05;
            c10r = 1:5;
            c01r = 1:5;
            range = Utils.createCombinations(0,c10r,c01r,pmaxr,lamr);
            idx = [5 9 10 13 14];
            prefix = 'zero_muscle_aniso_mooneytest';
            
            %% Test with no anisotropic muscle force and finer mr-coeffs
            % Did not work well at all (long comp times)
%             pmaxr = 400;
%             lamr = 2.05;
%             c10r = [0 logspace(-3,1,5)];
%             c01r = [0 logspace(-3,1,5)];
%             range = Utils.createCombinations(0,c10r,c01r,pmaxr,lamr);
%             idx = [5 9 10 13 14];
%             prefix = 'zero_muscle_aniso_mooneytest_withzero';

            %% Default run
%             range = double.empty(0,1);
%             idx = [];
%             prefix = 'default_params';
            
            %% -- EACH to be combinable with --
            
            %% Geoconfig 1: With activation
            % Straight fibres in x direction, "loose" ends that
            % allow the geometry to expand when compressed
            c = experiments.IsometricActivation('Tag',prefix,'BC',1,'FL',1,'GeoNr',Geo);
            cap = 'Movable setup with x-aligned fibres';
            
            %% Geoconfig 2
            % Straight fibres in x direction, but completely
            % fixed ends to ensure comparability with version 3
%             c = experiments.IsometricActivation('Tag',[prefix '_fixed'],'BC',2,'FL',1,'GeoNr',Geo);
%             cap = 'Fixed setup with x-aligned fibres';
             
            %% Geoconfig 3
            % Diagonal fibres in xz direction with completely
            % fixed ends
%             c = experiments.IsometricActivation('Tag',[prefix '_fixed_xzfibre'],'BC',1,'FL',1,'GeoNr',Geo);
%             cap = 'Fixed setup with xz-diagonal fibres';

            %% Geoconfig 4: No activation (fit for mooney-rivlin)
            % Straight fibres in x direction, "loose" ends that
            % allow the geometry to expand when compressed
%             c = experiments.IsometricActivation('Tag',prefix,'BC',1,'FL',1,'GeoNr',Geo,'Activate',0);
%             cap = 'Movable setup with x-aligned fibres, no activation';
            
            %% NOT CONFIGURABLE PART
            m = muscle.Model(c);
            e = tools.ExperimentRunner(m);
            mus = repmat(m.DefaultMu,1,size(range,2));
            mus(idx,:) = range;
            data = e.runExperimentsCached(mus);
            
            sp = c.ExperimentalStretchMillimeters;
            [sp,idx] = sort(sp);
            passive = data.o(:,idx,2);
            %passive(:,sp < 0) = -passive(:,sp < 0);
            active = data.o(:,idx,1);%-passive;
            
            pm = PlotManager;
            pm.LeaveOpen = true;
            pm.UseFileTypeFolders = false;
            tv = c.TargetOutputValues(idx,:);
            nex = size(mus,2);
            % Select best mooney-rivlin parameter for muscle
            if ~c.Options.Activate
                orig = repmat(tv(:,2),1,nex);
                diff = Norm.L2(-passive' - orig);
%                 reldiff = Norm.L2((-passive' - orig)./abs(orig));
%                 reldiff = sort(mean(abs((-passive' - orig)./abs(orig))));
                ax = pm.nextPlot(['fit_passive_mooney_rivlin_' c.Options.Tag],...
                    'Absolute L2 fitting errors','stretch','error');
                %plot(ax,sp,tv(:,2),'rx');
                %hold(ax,'on');
                plot(ax,1:nex,diff,'rx');
                [~, sel] = sort(diff);
                sel = sel(1:3);
            else
                sel = 1:nex;
            end
            
            %% Normal output plot
            ax = pm.nextPlot(['isomet_' c.Options.Tag],...
                ['Experiment: ' cap],'stretch','forces [mN]');
            plot(ax,sp,-passive(sel,:)','b');
            hold(ax,'on');
            plot(ax,sp,tv(:,2),'bx');
            if c.Options.Activate
                plot(ax,sp,-active');
                plot(ax,sp,tv(:,1)-tv(:,2),'rx');
            end
            
            % Set label
            % compute percent
            w = c.PosFE.Geometry.Width;
            sp_perc = (((w+sp)./w)-1)*100;
            lab = sprintf('%dmm (%.3g%%)|',reshape([sp; sp_perc],1,[]));
            set(ax,'XTick',sp,'XTickLabel', lab(1:end-1));
            
            disp('Values:');
            disp(active);
            disp('Ratios:');
            disp(bsxfun(@rdivide,active,active(1,:)));
            disp('Computation times [s]:');
            disp(data.ct);
            
            pm.savePlots(c.ImgDir,'Format',{'eps','jpg'});
        end
    end
    
end

