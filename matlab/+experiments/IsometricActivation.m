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
    
    methods
        function this = IsometricActivation(varargin)
            this = this@experiments.AExperimentModelConfig(varargin{:});
            this.addOption('BC',1);
            this.init;
            
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
            
            m.DefaultMu(5) = .001;
            m.DefaultMu(6) = 40;
            m.DefaultMu(7) = 1;
            m.DefaultMu(8) = 40;
            
            m.DefaultMu(13) = 380; % [kPa]
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
        
        function tmr = getTendonMuscleRatio(this, points)
            % Returns the [0,1] ratio between tendon and muscle at all
            % specified points
            %
            % This method simply returns an all-zero ratio, meaning muscle only. 
            tmr = zeros(1,size(points,2));
            o = this.Options;
            if o.GeoNr == 2
                tmr = this.TMRFun(points(2,:),points(3,:));
            end            
        end
        
        function tmr = TMRFun(~, y, z)
            tmr = zeros(size(y));
            fr = @(x)1.2*sqrt(max(0,27-x))-2;
            right = z > fr(y);
            tmr(right) = 1;
            fl = @(x)-1*sqrt(x)+1;
            left = z < fl(y);
            tmr(left) = 1;
        end
        
        function o = getOutputOfInterest(this, t, y)
            m = this.Model;
            % Get the residual dirichlet forces
            df = m.getResidualForces(t,y);
            
            % Get the position at which to determine the passive forces
            passive_pos = floor((this.PositioningTime + this.RelaxTime)/m.dt);
            % Get the range at which to determine the max forces
            act_range = passive_pos + 1 : size(y,2);
            switch this.Options.GeoNr
                case 1
                    idx = m.getVelocityDirichletBCFaceIdx(3,2);        
                case 2    
                    %this.Geometry.Nodes(2,:) > 33 & m.System.bool_v_bc_nodes
                    idx = m.getVelocityDirichletBCFaceIdx(37:48,4);
            end
            o(1) = max(abs(sum(df(idx,act_range),1)));
            o(2) = sum(df(idx,passive_pos));            
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
            if this.Options.GeoNr == 1
                [pts, cubes] = geometry.Cube27Node.DemoGrid(linspace(0,25,4),[0 7],[0 4]);
                geo = geometry.Cube27Node(pts,cubes);
            else
                cent = 12;
                k = kernels.GaussKernel(10);
                k2 = kernels.GaussKernel(12);
                rad = @(t)[k.evaluate(t,cent)*3.5 k2.evaluate(t,cent)*2]';
                geo = Belly.getBelly(4,35,'Radius',rad,'Layers',[.8 1],'InnerRadius',.3);
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
                case 2
                    displ_dir(2,geo.Nodes(2,:) == 0) = true;
                    displ_dir(:,1) = true;
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
                case 2
                    velo_dir(2,geo.Nodes(2,:) > 33) = true;
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
                case 2
                    anull(2,:,:) = 1;
                    anull(3,:,:) = 1;
            end
        end
    end
    
    methods(Static)
        
        function runExperiments()
            % Runs the series of isometric tests for different Pmax values.
            % Several variants can be chosen:
            %
            
            %% Mus 1 - Only PMAX
            range = 320:20:400;
            idx = 13;
            mustr = sprintf('-%d',range);
            prefix = 'pmax';
            
            %% Mus 1 - Only lambda_opt
%             range = 2:.1:2.2;
%             idx = 14;
%             mustr = sprintf('-%g',range);
%             prefix = 'lamopt';
            
            %% Mus 3 - pmax/lamdaopt
%             range = Utils.createCombinations(400:10:460,2:.025:2.2);
%             idx = [13 14];
%             mustr = [sprintf('-%g',400:10:460) '/' sprintf('-%g',2:.025:2.2)];
%             prefix = 'pmax_lamopt';
            
            %% -- EACH to be combinable with --
            
            %% Geoconfig 1
            % Straight fibres in x direction, "loose" ends that
            % allow the geometry to expand when compressed
            c = experiments.IsometricActivation('Tag',prefix,'BC',1,'FL',1);
            cap = 'Movable setup with x-aligned fibres';
            
            %% Geoconfig 2
            % Straight fibres in x direction, but completely
            % fixed ends to ensure comparability with version 3
%             c = experiments.IsometricActivation('Tag',[prefix '_fixed'],'BC',2,'FL',1);
%             cap = 'Fixed setup with x-aligned fibres';
             
            %% Geoconfig 3
            % Diagonal fibres in xz direction with completely
            % fixed ends
%             c = experiments.IsometricActivation('Tag',[prefix '_fixed_xzfibre'],'BC',3,'FL',1);
%             cap = 'Fixed setup with xz-diagonal fibres';
            
            %% NOT CONFIGURABLE PART
            m = muscle.Model(c);
            e = tools.ExperimentRunner(m);
            mus = repmat(m.DefaultMu,1,size(range,2));
            mus(idx,:) = range;
            data = e.runExperimentsCached(mus);
            
            sp = c.ExperimentalStretchMillimeters;
            [sp,idx] = sort(sp);
            passive = data.o(:,idx,2);
            passive(:,1:5) = -passive(:,1:5);
            active = data.o(:,idx,1)-passive;
            
            pm = PlotManager;
            pm.LeaveOpen = true;
            pm.UseFileTypeFolders = false;
            ax = pm.nextPlot(['isomet_' c.Options.Tag],...
                ['Pmax' mustr ' experiment: ' cap],'stretch','forces [mN]');
            plot(ax,sp,active','r',sp,passive','b');
            hold(ax,'on');
            tv = c.TargetOutputValues(idx,:);
            plot(ax,sp,tv(:,1)-tv(:,2),'rx',sp,tv(:,2),'bx');
            
            % Set label
            % compute percent
            w = c.PosFE.Geometry.Width;
            sp_perc = (((w+sp)./w)-1)*100;
            lab = sprintf('%dmm (%g%%)|',reshape([sp; sp_perc],1,[]));
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

