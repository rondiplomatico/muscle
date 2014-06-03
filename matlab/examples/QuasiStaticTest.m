classdef QuasiStaticTest < muscle.AModelConfig
% A long geometry with 20% deviation from default cubic positions and
% complex fibre structure

    properties(Constant)
        OutputDir = fullfile(fileparts(which(mfilename)),'quasistatictest');
    end

    properties
        Case;
        GeoNr;
    end

    methods
        function this = QuasiStaticTest(thecase, geonr)
            if nargin < 2
                geonr = 1;
                if nargin < 1
                    thecase = 1;
                end
            end
            switch geonr
                case 1
                    [pts, cubes] = geometry.Cube8Node.DemoGrid(0:20:40,[0 20],[0 20]);
                     geo = geometry.Cube8Node(pts, cubes);
                case 2
                    geo = Belly.getBelly(8, 50, 3, .5, 15);
                case 3
                    s = load(fullfile(fileparts(which(mfilename)),'..','CMISS','EntireTA.mat'));
                    geo = s.geo27;
            end
            this = this@muscle.AModelConfig(geo);
            this.Case = thecase;
            this.GeoNr = geonr;
        end
        
        function configureModel(this, m)
            os = m.ODESolver;
            switch this.GeoNr
                case 1
                    os.RelTol = .01;
                    os.AbsTol = .05;
                case {2,3}
                    os.RelTol = .01;
                    os.AbsTol = .08;
            end
            m.EnableTrajectoryCaching = true;
        end
        
        function o = getOutputOfInterest(this, m, t, uvw)
            geo = m.Config.PosFE.Geometry;
            if this.Case == 1
                uvw = m.System.includeDirichletValues(t, uvw);
                switch this.GeoNr
                    case 1
                        facenode_idx = m.getFaceDofsGlobal(2,2,1);
                        o(1,:) = mean(uvw(facenode_idx,:),1);
                        o(2,:) = mean(uvw(facenode_idx+geo.NumNodes*3,:),1);
                    case 2
                        facenode_idx = [];
                        for k = 1:4
                            facenode_idx = [facenode_idx; m.getFaceDofsGlobal(k,3,2)];%#ok
                        end
                        % Save some work
                        facenode_idx = unique(facenode_idx);
                        o(1,:) = mean(uvw(facenode_idx,:),1);
                        o(2,:) = mean(uvw(facenode_idx+geo.NumNodes*3,:),1);
                    case 3
                        facenode_xidx = m.getFaceDofsGlobal(6,3,1);
                        facenode_yidx = m.getFaceDofsGlobal(6,3,2);
                        facenode_xidx = unique(facenode_xidx);
                        facenode_yidx = unique(facenode_yidx);
                        
                        % Store x,y positions and velocities separately
                        o(3,:) = mean(uvw(facenode_xidx,:),1);
                        o(4,:) = mean(uvw(facenode_yidx,:),1);
                        o(5,:) = mean(uvw(facenode_xidx+geo.NumNodes*3,:),1);
                        o(6,:) = mean(uvw(facenode_yidx+geo.NumNodes*3,:),1);
                        
                        % Return norm of difference vector (to initial pos)
                        o(1,:) = Norm.L2(o([3 4],:) - o([3 4],1)*ones(1,size(o,2)));
                        o(2,:) = Norm.L2(o([5 6],:) - o([5 6],1)*ones(1,size(o,2)));
                end
            elseif this.Case == 2
                df = m.getResidualForces(t, uvw);
                switch this.GeoNr
                    case 1
                        idx = m.getDirichletBCFaceIdx(2,2);
                        o(1,:) = mean(df(idx,:),1);
                    case 2
                        xidx = m.getDirichletBCFaceIdx(1:4,3,1);
                        yidx = m.getDirichletBCFaceIdx(1:4,3,2);
                        zidx = m.getDirichletBCFaceIdx(1:4,3,3);
                        vec = [mean(df(xidx,:),1); mean(df(yidx,:),1); mean(df(zidx,:),1)];
                        o(1,:) = Norm.L2(vec);
                        o(2:4,:) = vec;
                    case 3
                        xidx = m.getDirichletBCFaceIdx(6,3,1);
                        yidx = m.getDirichletBCFaceIdx(6,3,2);
                        zidx = m.getDirichletBCFaceIdx(6,3,3);
                        vec = [mean(df(xidx,:),1); mean(df(yidx,:),1); mean(df(zidx,:),1)];
                        o(1,:) = Norm.L2(vec);
                        o(2:4,:) = vec;
                end
            end
        end
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
                case 2
                    % Fix back side
                    for k = geo.NumElements-3:geo.NumElements
                        displ_dir(:,geo.Elements(k,geo.MasterFaces(4,:))) = true;
                    end
                case 3
                    % Fix broad end of TA
                    displ_dir(:,geo.Elements(8,geo.MasterFaces(4,:))) = true;
                    displ_dir(:,geo.Elements(8,geo.MasterFaces(2,:))) = true;
            end
            %% In case 2, we additionally fix the opposite side (and then activate)
            if this.Case == 2
               switch this.GeoNr
                   case 1
                       % Fix back
                       displ_dir(:,geo.Elements(2,geo.MasterFaces(2,:))) = true;
                   case 2
                       % Fix fron side
                       for k = 1:4
                           displ_dir(:,geo.Elements(k,geo.MasterFaces(3,:))) = true;
                       end
                   case 3
                       % Fix thin end of TA
                       displ_dir(:,geo.Elements(6,geo.MasterFaces(3,:))) = true;
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
        function runTestCase1(geonr)
            if nargin < 1
                geonr = 1;
            end
            mc = QuasiStaticTest(1,geonr);
            
            % The activation rates (alpha increase per ms)
            ramptimes = [.1:.1:1 1:.2:2 3:10 20 40 50 100 300 600 1000 2000 8000 60000]; % [ms]
            
%             pm = PlotManager(false,2,2);
            pm = PlotManager;
            pm.AutoTickMarks = 0;
            pm.LeaveOpen = true;
            
            c = ColorMapCreator;
            c.useJet([0.01 0.05 0.1 .5]);
%             c.addColor(0.005,[0 .7 0]);
%             c.addColor(0.01,[.7 .5 .2],[.5 .7 .2]);
%             c.addColor(0.05,[1 1 .4],[.7 .5 .2]*.5);
%             c.addColor(0.1,[.7 0 0],[1 1 .4]*.3);
%             c.addColor(10,[1 0 1]);
           
%             visc = [0.001 0.01 .1 1 10];
            visc = [0.01 .1 1 10];
            nvisc = length(visc);
            for k=1:nvisc
                m = muscle.Model(mc);
                v = visc(k);
                m.System.Viscosity = v;
                file = fullfile(QuasiStaticTest.OutputDir,sprintf('case1_geo%d_visc%g.mat',geonr,visc(k)));
                
                if exist(file,'file') == 2
                    load(file);
                else
                    [alphavals, pos] = QuasiStaticTest.runAlphaRamp(m, ramptimes); %globaltimes, globalpos
                    save(file,'alphavals','m','pos','ramptimes');
%                     save(file,'alphavals','m','pos','globaltimes','globalpos','ramptimes');
                end
                
                ec = [.3 .3 .3];
                h = pm.nextPlot(sprintf('pos_visc%g',v),...
                    sprintf('Averaged X-position [mm] of right face for different rates and activation level\nViscosity=%g',v),...
                    'Activation increase rate [1/ms]','alpha [-]');
                surf(h,1./ramptimes, alphavals, pos(:,:,1)','FaceColor','interp','EdgeColor',ec);
                set(h,'XScale','log');
                view(-150, 45);
                
                h = pm.nextPlot(sprintf('velo_visc%g',v),...
                    sprintf('Averaged X-velocity [mm/ms] of right face for different rates and activation level\nViscosity=%g',v),...
                    'Activation increase rate [1/ms]','alpha [-]');
                surf(h,1./ramptimes, alphavals, abs(pos(:,:,2))','FaceColor','interp','EdgeColor',ec);
                set(h,'XScale','log');
                view(-136, 56);

                h = pm.nextPlot(sprintf('pos_abserr_visc%g',v),...
                    sprintf('Absolute error in of X-position [mm] for different rates and activation level\nagainst quasi-static positions (rate=%g/ms)\nViscosity=%g',1/ramptimes(end),v),...
                    'Activation increase rate [1/ms]','alpha [-]');
                staticpos = repmat(pos(end,:,1),size(pos,1),1);
                diffX = abs(pos(:,:,1) - staticpos)';
                surf(h,1./ramptimes, alphavals, diffX,'FaceColor','interp','EdgeColor',ec);
                set(h,'XScale','log');
                view(-150, 45);
                
                h = pm.nextPlot(sprintf('pos_abserr_topview_visc%g',v),...
                    sprintf('Absolute error of X-position [mm] for different rates and activation level\nagainst quasi-static positions (rate=%g/ms)\nViscosity=%g',1/ramptimes(end),v),...
                    'Activation increase rate [1/ms]','alpha [-]');
                surf(h,1./ramptimes, alphavals, diffX,'FaceColor','interp','EdgeColor',ec);
                set(h,'XScale','log'); colorbar;
                view(h,-180, 90);
                
                h = pm.nextPlot(sprintf('pos_relerr_visc%g',v),...
                    sprintf('Relative error of X-position [mm] for different rates and activation level\nagainst quasi-static positions (rate=%g/ms)\nViscosity=%g',1/ramptimes(end),v),...
                    'Activation increase rate [1/ms]','alpha [-]');
                Z = abs(diffX ./ staticpos');
                surf(h,1./ramptimes, alphavals, Z ,'FaceColor','interp','EdgeColor',ec);
                c.LogPlot = false;
                colormap(h,c.create(Z)); colorbar;
                set(h,'XScale','log');
                view(-150, 45);

                nsteps = length(alphavals)-1;
                dt = meshgrid(ramptimes, 1:nsteps)';
                velo = diff(staticpos,[],2) ./ (dt/nsteps);
                h = pm.nextPlot(sprintf('velo_estim_visc%g',v),...
                    sprintf('Inferred velocities [mm/ms] for different rates and "quasi time step"\nViscosity=%g',v),...
                    '\Delta t [1/ms]','alpha  * \Delta t ');
                Z = abs(velo)';
                LogPlot.logsurfc(h,1./ramptimes, alphavals(2:end), Z,'FaceColor','interp','EdgeColor',ec);
                set(h,'XScale','log'); axis(h,'tight');
                view(-150, 45);
                
                h = pm.nextPlot(sprintf('velo_relerr_visc%g',v),...
                    sprintf('Relative error between computed (via rate) velocities [mm/ms] for different rates and "quasi time step"\nViscosity=%g',v),...
                    '\Delta t [1/ms]','Fraction of timestep');
                Z = abs((velo-pos(:,2:end,2))./pos(:,2:end,2))';
                LogPlot.logsurfc(h,1./ramptimes, alphavals(2:end), Z,'FaceColor','interp','EdgeColor',ec);
                set(h,'XScale','log'); axis(h,'tight');
                c.LogPlot = true;
                colormap(h,c.create(Z)); colorbar;
                view(-125, 50);
                
                h = pm.nextPlot(sprintf('velo_relerr_maxmean_visc%g',v),...
                    sprintf('Relative error between computed (via rate) velocities [mm/ms] for different rates\nMax/Mean over all "quasi time steps"\nViscosity=%g',v),...
                    '\Delta t [1/ms]','Error');
                maxv = max(abs((velo-pos(:,2:end,2))./pos(:,2:end,2)),[],2);
                meanv = mean(abs((velo-pos(:,2:end,2))./pos(:,2:end,2)),2);
                loglog(h,1./ramptimes, [maxv meanv]); axis(h,'tight');
                legend('Max relative error','Mean relative error','Location','NorthWest');
            end
            
            m.plotGeometrySetup(pm);
            
            pm.done;
            pm.FilePrefix = sprintf('case1_geo%d',geonr);
            pm.ExportDPI = 200;
            pm.SaveFormats = {'jpg'};
            pm.savePlots(QuasiStaticTest.OutputDir,'Close',true);
        end
        
        function runTestCase2(geonr)
            % The 
            if nargin < 1
                geonr = 1;
            end
            mc = QuasiStaticTest(2,geonr);
            
            % The activation rates (alpha increase per ms)
            ramptimes = [.1:.1:1 1:.2:2 3:10 20 40 50 100 300 1000 10000 60000]; % [ms]
            
%             pm = PlotManager(false,2,2);
            pm = PlotManager;
            pm.AutoTickMarks = 0;
            pm.LeaveOpen = true;
            c = ColorMapCreator;
            c.useJet([0.1 0.2 .5 1]);
            
            visc = [0.001 0.01 .1 1 10];
%             visc = [0.01 .1 1 10];
            nvisc = length(visc);
            for k=1:nvisc
                v = visc(k);
                file = fullfile(QuasiStaticTest.OutputDir,sprintf('case2_geo%d_visc%g.mat',geonr,visc(k)));
                if exist(file,'file') == 2
                    load(file);
%                     [alphavals, forces] = QuasiStaticTest.runAlphaRamp(m, ramptimes);
%                     save(file,'forces','-APPEND');
                else
                    m = muscle.Model(mc);
                    m.System.Viscosity = v;
                    [alphavals, forces] = QuasiStaticTest.runAlphaRamp(m, ramptimes);
                    save(file,'alphavals','m','forces','ramptimes');
                end
                
                ec = [.3 .3 .3];
                h = pm.nextPlot(sprintf('forcenorm_visc%g',v),...
                    sprintf('Force [N] (averaged over one fixed side)\nViscosity=%g',v),...
                    'Activation increase rate [1/ms]','alpha [-]');
                surf(h,1./ramptimes, alphavals, forces(:,:,1)','FaceColor','interp','EdgeColor',ec);
                set(h,'XScale','log');
                view(-150, 45);
                
                staticnorm = repmat(forces(end,:,1),size(forces,1),1);
                diff = staticnorm - forces(:,:,1);
                h = pm.nextPlot(sprintf('forcenorm_visc%g',v),...
                    sprintf('Error of Force [N] compared to quasi-static solution (rate %g/ms)\nViscosity=%g',1/ramptimes(end),v),...
                    'Activation increase rate [1/ms]','alpha [-]');
                surf(h,1./ramptimes, alphavals, diff','FaceColor','interp','EdgeColor',ec);
                set(h,'XScale','log');
                view(-150, 45);
                
                h = pm.nextPlot(sprintf('forcenorm_visc%g',v),...
                    sprintf('Relative error of Force [N] compared to quasi-static solution (rate %g/ms)\nViscosity=%g',1/ramptimes(end),v),...
                    'Activation increase rate [1/ms]','alpha [-]');
                Z = abs(diff ./ staticnorm)';
                surf(h,1./ramptimes, alphavals, Z,'FaceColor','interp','EdgeColor',ec);
                set(h,'XScale','log');
%                 c.LogPlot = true;
                colormap(h,c.create(Z)); colorbar;
                view(-150, 45);
                
                %% Single dimensions
%                 h = pm.nextPlot(sprintf('xpos_rate_alpha_visc%g',v),...
%                     sprintf('X\nViscosity=%g',v),...
%                     'Activation increase rate [1/ms]','alpha [-]');
%                 surf(h,1./ramptimes, alphavals, forces(:,:,2)','FaceColor','interp','EdgeColor',ec);
%                 set(h,'XScale','log');
%                 view(-150, 45);
%                 
%                 h = pm.nextPlot(sprintf('xpos_rate_alpha_visc%g',v),...
%                     sprintf('Y\nViscosity=%g',v),...
%                     'Activation increase rate [1/ms]','alpha [-]');
%                 surf(h,1./ramptimes, alphavals, forces(:,:,3)','FaceColor','interp','EdgeColor',ec);
%                 set(h,'XScale','log');
%                 view(-150, 45);
%                 
%                 h = pm.nextPlot(sprintf('xpos_rate_alpha_visc%g',v),...
%                     sprintf('Z\nViscosity=%g',v),...
%                     'Activation increase rate [1/ms]','alpha [-]');
%                 surf(h,1./ramptimes, alphavals, forces(:,:,4)','FaceColor','interp','EdgeColor',ec);
%                 set(h,'XScale','log');
%                 view(-150, 45);

            end
            
            m.plotGeometrySetup(pm);
            
%             pm.done;
%             pm.FilePrefix = sprintf('case2_geo%d',geonr);
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

