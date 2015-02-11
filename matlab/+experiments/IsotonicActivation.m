classdef IsotonicActivation < experiments.AExperimentModelConfig
    % Implements the isometric contraction experiment
    
    properties(Constant)
        ExperimentalStretchMillimeters = [0 1 -1 2 -2 3 -4 4 -5 5 -6];
    end
    
    properties
        ActivationTime = 100; %[ms]
        PositioningTime = 50; %[ms]
        RelaxTime = 10; %[ms]
        
        PreStretchLength = 10; %[mm]
    end
    
    methods
        function this = IsotonicActivation(varargin)
            this = this@experiments.AExperimentModelConfig(varargin{:});
            %this.addOption('BC',1);
            this.init;
            
            % We need computed initial conditions on this one
            this.RequiresComputedInitialConditions = true;
            
            this.NumOutputs = 2;
            this.NumConfigurations = 1;
            this.TargetOutputValues = 1;
            
            this.VelocityBCTimeFun = tools.ConstantUntil(this.PositioningTime,.01);
        end
        
        function configureModel(this, m)
            configureModel@muscle.AModelConfig(this, m);
            os = m.ODESolver;
            os.RelTol = 1e-4;
            os.AbsTol = .001;
            
            m.T = this.PositioningTime + this.RelaxTime + this.ActivationTime;
            m.dt = m.T / 300;
            
            m.DefaultMu(5) = 1;
            m.DefaultMu(6) = 20;
            m.DefaultMu(7) = 1;
            m.DefaultMu(8) = 25;
            
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
        
        function o = getOutputOfInterest(this, t, y)
            m = this.Model;
            % Get the residual dirichlet forces
            df = m.getResidualForces(t,y);
            
%             m.plot(t,y,'DF',df,'F',4);
            % Get the position at which to determine the passive forces
            passive_pos = floor((this.PositioningTime + this.RelaxTime)/m.dt);
            % Get the range at which to determine the max forces
            act_range = passive_pos + 1 : size(y,2);
            switch this.Options.GeoNr
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
        
        function geo = getGeometry(this)
            if this.Options.GeoNr == 1
                [pts, cubes] = geometry.Cube27Node.DemoGrid(linspace(0,25,4),[0 7],[0 4]);
                geo = geometry.Cube27Node(pts,cubes);
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
                    velo_dir_val(velo_dir) = velo;
                case 2
            end
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
            end
        end
    end
    
    methods(Static)
        
        function runExperiments()
            % Runs the series of isometric tests for different Pmax values.
            % Several variants can be chosen:
            %
            
            %% Mus 1 - Only PMAX
%             range = 320:20:400;
%             idx = 13;
%             mustr = sprintf('-%d',range);
%             prefix = 'pmax';
            
            %% Mus 1 - Only lambda_opt
%             range = 2:.1:2.2;
%             idx = 14;
%             mustr = sprintf('-%g',range);
%             prefix = 'lamopt';
            
            %% Mus 3 - pmax/lamdaopt
            range = Utils.createCombinations(400:10:460,2:.025:2.2);
            idx = [13 14];
            mustr = [sprintf('-%g',400:10:460) '/' sprintf('-%g',1:.05:1.2)];
            prefix = 'pmax_lamopt';
            
            %% -- EACH to be combinable with --
            
            %% Geoconfig 1
            % Straight fibres in x direction, "loose" ends that
            % allow the geometry to expand when compressed
%             c = experiments.IsometricActivation('Tag',prefix,'BC',1,'FL',1);
%             cap = 'Movable setup with x-aligned fibres';
            
            %% Geoconfig 2
            % Straight fibres in x direction, but completely
            % fixed ends to ensure comparability with version 3
%             c = experiments.IsometricActivation('Tag',[prefix '_fixed'],'BC',2,'FL',1);
%             cap = 'Fixed setup with x-aligned fibres';
             
            %% Geoconfig 3
            % Diagonal fibres in xz direction with completely
            % fixed ends
            c = experiments.IsometricActivation('Tag',[prefix '_fixed_xzfibre'],'BC',3,'FL',1);
            cap = 'Fixed setup with xz-diagonal fibres';
            
            %% NOT CONFIGURABLE PART
            m = muscle.Model(c);
            e = experiments.ExperimentRunner(m);
            mus = repmat(m.DefaultMu,1,size(range,2));
            mus(idx,:) = range;
            data = e.runExperimentsCached(mus);
            
            sp = c.ExperimentalStretchMillimeters;
            [sp,idx] = sort(sp);
            passive = data.o(:,idx,2);
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

