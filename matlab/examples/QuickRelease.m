classdef QuickRelease < muscle.AModelConfig
% A long geometry with 20% deviation from default cubic positions and
% complex fibre structure

    properties(Constant)
        OutputDir = fullfile(fileparts(which(mfilename)),'quickrelease');
    end

    properties
        GeoNr;
        ICCompMode;
    end
    
    properties(SetAccess=private)
        % The time the velocity BCs are applied in order to reach the
        % initial condition
        icMovetime = 200; % [ms]
        
        % The time used to increase the muscle activation
        icAlphaRampTime = 100; % [ms]
        
        % The time after movement and activation is fully done until the
        % ICs are extracted ("quasi static")
        icRelaxTime = 300; % [ms]
    end

    methods
        function this = QuickRelease(geonr, iccomp)
            Utils.ensureDir(QuickRelease.OutputDir);
            if nargin < 2
                iccomp = false;
                if nargin < 1
                    geonr = 2;
                end
            end
            switch geonr
                case 1
                    [pts, cubes] = geometry.Cube8Node.DemoGrid(0:20:40,[0 20],[0 20]);
                    geo = geometry.Cube8Node(pts, cubes);
                case 2
                    geo = Belly.getBelly(3, 50, 3, .5, 15);
                case 3
                    s = load(fullfile(fileparts(which(mfilename)),'..','CMISS','EntireTA.mat'));
                    geo = s.geo27;
            end
            this = this@muscle.AModelConfig(geo);
            this.ICCompMode = iccomp;
            this.GeoNr = geonr;
        end
        
        function configureModel(this, m)
            os = m.ODESolver;
            switch this.GeoNr
                case 1
                    if this.ICCompMode
                        os.RelTol = .001;
                        os.AbsTol = .01;
                    else
                        os.RelTol = .1;
                        os.AbsTol = .1;
                    end
                case {2,3}
                    os.RelTol = .01;
                    os.AbsTol = .08;
            end
            % For IC comps, use an alpha ramp. Otherwise, start with full
            % activation straight away
            f = m.System.f;
            alpha = .1;
            if this.ICCompMode
                f.alpha = this.getAlphaRamp(this.icAlphaRampTime, alpha);
                m.T = max(this.icMovetime, this.icAlphaRampTime) + this.icRelaxTime;
                m.dt = m.T / 100;
                m.System.ApplyVelocityBCUntil = this.icMovetime;
            else
                m.T = 100;
                m.dt = .1;
                f.alpha = @(t)alpha;
                m.EnableTrajectoryCaching = true;
            end
        end
        
        function x0 = getX0(this, x0)
            if ~this.ICCompMode
                s = load(fullfile(QuickRelease.OutputDir,sprintf('ic_geo%d.mat',this.GeoNr)));
                x0 = s.x0;
            end
        end
        
%         function o = getOutputOfInterest(this, m, t, uvw)
%             geo = m.Config.PosFE.Geometry;
%             
%            
%         end
        
        
    end
    
    methods(Access=protected)
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            
            geo = this.PosFE.Geometry;
            %% Dirichlet conditions: Position (fix one side)
            % This is done for each test case.
            switch this.GeoNr
                case 1
                    % Fix front
                    displ_dir(:,geo.Elements(1,geo.MasterFaces(1,:))) = true;
%                 case 2
%                     % Fix back side
%                     for k = geo.NumElements-3:geo.NumElements
%                         displ_dir(:,geo.Elements(k,geo.MasterFaces(4,:))) = true;
%                     end
%                 case 3
%                     % Fix broad end of TA
%                     displ_dir(:,geo.Elements(8,geo.MasterFaces(4,:))) = true;
%                     displ_dir(:,geo.Elements(8,geo.MasterFaces(2,:))) = true;
            end
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(this, velo_dir, velo_dir_val)
            % Determines the dirichlet velocities.
            %
            % The unit for the applied quantities is [mm/ms] = [m/s]
            %
            % These conditions are required to compute the initial
            % conditions for the quick release test.
            if this.ICCompMode
                geo = this.PosFE.Geometry;
                box = geo.getBoundingBox;
                f = this.Model.System.f;
                switch this.GeoNr
                    case 1
                        xdia = box(2)-box(1);
                        %totaldistance = (f.lambdafopt-1)*xdia;
                        totaldistance = .01*xdia;
                        velo_dir(1,geo.Elements(2,geo.MasterFaces(2,:))) = true;
                        velo_dir_val(velo_dir) = totaldistance/this.icMovetime;
                end
            end
        end
        
        function anull = seta0(this, anull)
            switch this.GeoNr
                case {1, 3}
                    % Fibres in x direction
                    anull(1,:,:) = 1;
                case 2
                    % Fibres in y direction
                    anull(2,:,:) = 1;
            end
        end
    end
    
    methods(Static)
        function runTest(geonr)
            if nargin < 1
                geonr = 1;
            end
            
            %% Initial conditions
            file = fullfile(QuickRelease.OutputDir,sprintf('ic_geo%d.mat',geonr));
            if exist(file,'file') ~= 2
                mc = QuickRelease(geonr, true);
                m = muscle.Model(mc);
                [t, y] = m.simulate(0);
                m.plot(t,y);
                y = m.System.includeDirichletValues(t,y);
                x0 = y(:,end);%#ok
                save(file,'x0');
            end
            
            mc = QuickRelease(geonr, false);
            
            pm = PlotManager(false,2,2);
%             pm = PlotManager;
            pm.AutoTickMarks = 0;
            pm.LeaveOpen = true;
            
            c = ColorMapCreator;
            c.useJet([0.01 0.05 0.1 .5]);
            
            visc = 10;
%             visc = [0.001 0.01 .1 1 10];
            nvisc = length(visc);
            for k=1:nvisc
                m = muscle.Model(mc);
                v = visc(k);
                m.System.Viscosity = v;
                file = fullfile(QuickRelease.OutputDir,sprintf('geo%d_visc%g.mat',geonr,visc(k)));
                
                if exist(file,'file') == 2
                    load(file);
                else
                    [t, y, ct] = m.simulate(0);
                    m.plot(t,y);
%                     o = m.Config.getOutputOfInterest(m, t, y);
%                     save(file,'m','y','o','ct');
                end
                
%                 ec = [.3 .3 .3];
%                 h = pm.nextPlot(sprintf('pos_visc%g',v),...
%                     sprintf('Averaged X-position [mm] of right face for different rates and activation level\nViscosity=%g',v),...
%                     'Activation increase rate [1/ms]','alpha [-]');
%                 surf(h,1./ramptimes, alphavals, pos(:,:,1)','FaceColor','interp','EdgeColor',ec);
%                 set(h,'XScale','log');
%                 view(-150, 45);
%                 
%                 h = pm.nextPlot(sprintf('velo_visc%g',v),...
%                     sprintf('Averaged X-velocity [mm/ms] of right face for different rates and activation level\nViscosity=%g',v),...
%                     'Activation increase rate [1/ms]','alpha [-]');
%                 surf(h,1./ramptimes, alphavals, abs(pos(:,:,2))','FaceColor','interp','EdgeColor',ec);
%                 set(h,'XScale','log');
%                 view(-136, 56);
            end
            
%             m.plotGeometrySetup(pm);
%             
%             pm.done;
%             pm.FilePrefix = sprintf('case1_geo%d',geonr);
%             pm.ExportDPI = 200;
%             pm.SaveFormats = {'jpg'};
%             pm.savePlots(QuasiStaticTest.OutputDir,'Close',true);
        end
        
     
        
        function [alphavals, output] = runAlphaRamp(m, ramptimes) %globaltimes, globalout
            f = m.System.f;
            rampsteps = 30;
            alphamax = 1;
            alphavals = linspace(0,alphamax,rampsteps+1);
            ramppos = 1:2:rampsteps*2+1;
%             relaxperc = 1.3;
            nrates = length(ramptimes);

%             globaloutsampling = round(ramptimes(end)/4);
%             globaltimes = linspace(0,max(ramptimes)+relaxtime,globaloutsampling);
%             globalout = zeros(nrates,globaloutsampling,2);

            output = zeros(nrates, rampsteps+1, 2);
            pi = ProcessIndicator('Computing for %d different alpha increase rates',nrates,false,nrates);
            for tidx = 1:nrates
                ramptime = ramptimes(tidx);
                m.T = ramptime;%*relaxperc;
                m.dt = ramptime/rampsteps/2;

                f.alpha = @(t)alphamax * ((t>ramptime) + (t<=ramptime).*t/ramptime);

                [t,y] = m.simulate(0);
                o = m.Config.getOutputOfInterest(m, t, y);
                % Extract position on global time grid (plotting only)
%                 [~, idx] = min(abs(globaltimes-t(end))); idx = idx(1);
                for k = 1:size(o,1)
                    output(tidx,:,k) = o(k,ramppos);
%                     globalout(tidx,1:idx,k) = interp1(t,o(k,:),globaltimes(1:idx),'cubic');
%                     globalout(tidx,idx+1:end,k) = globalout(tidx,idx,k);
                end

                % Check correct alpha values
                if norm(f.alpha(t(ramppos)) - alphavals) > 1e-15
                    error('Boo.');
                end
                pi.step;
            end
            pi.stop;
        end
    end
    
end

