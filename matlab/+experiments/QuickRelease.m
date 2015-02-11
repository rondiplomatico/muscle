classdef QuickRelease < experiments.AExperimentModelConfig
% Provides a model config and test scripts for a quick release test.
%
% The material parameters used are suggested by Thomas Heidlauf, see his
% thesis.
  
    properties(SetAccess=private)
        % The time the velocity BCs are applied in order to reach the
        % initial condition
        icMovetime = 100; % [ms]
        
        % The time used to increase the muscle activation
        icAlphaRampTime = 100; % [ms]
        
        % The time after movement and activation is fully done until the
        % ICs are extracted ("quasi static")
        icRelaxTime = 100000; % [ms]
        
        Loads;
        Pressures;
    end

    methods
        function this = QuickRelease(varargin)
            this = this@experiments.AExperimentModelConfig(varargin{:});
            this.init;
        end
        
        function configureModel(this, m)
            os = m.ODESolver;
            switch this.GeoNr
                case 1
                    if this.ICCompMode
                        os.RelTol = 1e-7;
                        os.AbsTol = .009;
                    else
                        os.RelTol = .01;
                        os.AbsTol = .1;
                    end
                case 2
                    if this.ICCompMode
                        os.RelTol = 1e-7;
                        os.AbsTol = .009;
                    else
                        os.RelTol = .01;
                        os.AbsTol = .1;
                    end
                case 3
                    if this.ICCompMode
                        os.RelTol = 1e-7;
                        os.AbsTol = .02;
                    else
                        os.RelTol = .01;
                        os.AbsTol = .8;
                    end
            end
            %% Material setup
            % Material set (see main comment)
            m.DefaultMu(5) = 0.00355439810963035; % b1 [kPa]
            m.DefaultMu(6) = 12.660539325481963; % d1 [-]
            m.DefaultMu(9) = 6.352e-10; % c10 [kPa]
            m.DefaultMu(10) = 3.627; % c01 [kPa]
            m.DefaultMu(13) = 250; % [kPa]
            m.DefaultMu(14) = 1.2; % [-]
            
            if this.ICCompMode
                % No activation needed for position/IC computing
                m.DefaultMu(2) = 0;
%                 f.alpha = this.getAlphaRamp(this.icAlphaRampTime, alpha);
                m.T = max(this.icMovetime, this.icAlphaRampTime) + this.icRelaxTime;
                m.dt = m.T / 300;
                this.VelocityBCTimeFun = tools.ConstantUntil(this.icMovetime);
            else
                m.T = 20;
                m.dt = .1;
                m.DefaultMu(2) = 2*m.dt;
                m.EnableTrajectoryCaching = true;
            end
        end
        
        function x0 = getX0(this, x0)
            if ~this.ICCompMode
                s = load(fullfile(QuickRelease.OutputDir,sprintf('geo%d_ic.mat',this.GeoNr)));
                x0 = s.x0;
            end
        end
        
        function o = getOutputOfInterest(this, m, t, uvw)
            geo = m.Config.PosFE.Geometry;
            [df,nf] = m.getResidualForces(t,uvw);
            uvw = m.System.includeDirichletValues(t, uvw);
            switch this.GeoNr
                case 1
                    % Get average velocity of loose end
                    facenode_idx = m.getFaceDofsGlobal(2,2,1);
                    o(1,:) = mean(uvw(facenode_idx+geo.NumNodes*3,:),1);
                    
                    % Get force on that end [N]
                    o(2,:) = abs(sum(nf,1))/1000;
                    % Also store the dirichlet forces (maybe we'll need it)
                    % [N]
                    o(3,:) = sum(df,1)/1000;
                case 2
                    facenode_idx = [];
                    for k = 1:4
                        facenode_idx = [facenode_idx; m.getFaceDofsGlobal(k,3,2)];%#ok
                    end
                    % Save some work
                    facenode_idx = unique(facenode_idx);
                    o(1,:) = mean(uvw(facenode_idx+geo.NumNodes*3,:),1);
                    
                    % Get force on that end [N]
                    o(2,:) = abs(sum(nf,1))/1000;
                    % Also store the dirichlet forces (maybe we'll need it)
                    % [N]
                    o(3,:) = sum(df,1)/1000;
                case 3
                    % Get average velocity of loose end
                    facenode_idx = m.getFaceDofsGlobal(6,3,2);
                    o(1,:) = mean(uvw(facenode_idx+geo.NumNodes*3,:),1);
                    
                    % Get force on that end [N]
                    o(2,:) = abs(sum(nf,1))/1000;
                    % Also store the dirichlet forces (maybe we'll need it)
                    % [N]
                    o(3,:) = sum(df,1)/1000;
            end
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is kiloPascal [kPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            P = [];
            if ~this.ICCompMode
                switch this.GeoNr
                    case 1
                        if elemidx == 2 && faceidx == 2
                            P = 1;
                        end
                    case 2
                        if elemidx <= 4 && faceidx == 3
                            P = 1;
                        end
                    case 3
                        if elemidx == 6 && faceidx == 3
                            P = 1;
                        end
                end
            end
        end
        
        function u = getInputs(this)
            % loads in [g]
            loads = [0 100 2000];

            % convert to pressure:
            % [g]/1000 = [kg]
            % [kg]*[m/s²] = [N]
            % [N]*1000 = [mN]
            % [mN]/[mm²] = [kPa]
            switch this.GeoNr
                case 1
                    a = this.PosFE.getFaceArea(2,2);
                case 2
                    a = this.PosFE.getFaceArea(4,3);
                case 3
                    a = this.PosFE.getFaceArea(6,3);
            end
            m = this.Model;
            pressures = (loads/1000*m.Gravity)*1000/a; % [kPa]
            u = {};
            for lidx = 1:length(loads)
                pressure = pressures(lidx);
                u{lidx} = this.getAlphaRamp(2*m.dt,pressure);%#ok
            end
            this.Loads = loads;
            this.Pressures = pressures;
        end
        
    end
    
    methods(Access=protected)
        
        function geo = getGeometry(this)
            switch this.Options.GeoNr
                case 1
                    [pts, cubes] = geometry.Cube8Node.DemoGrid(0:20:40,[0 20],[0 20]);
                    geo = geometry.Cube8Node(pts, cubes);
                case 2
                    geo = Belly.getBelly(3, 50, 'Radius', 4.5, 'InnerRadius', 1.5, 'Gamma', 15);
                case 3
                    s = load(fullfile(fileparts(which(mfilename)),'..','CMISS','EntireTA.mat'));
                    geo = s.geo27;
            end
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            
            geo = this.PosFE.Geometry;
            %% Dirichlet conditions: Position (fix one side)
            % This is done for each test case.
            switch this.GeoNr
                case 1
                    % Fix front
                    face = geo.MasterFaces(1,:);
                    displ_dir(:,geo.Elements(1,face(5))) = true;
                    displ_dir(1,geo.Elements(1,face([1:4 6:9]))) = true;
                case 2
                    % Fix back side
                    for k = geo.NumElements-3:geo.NumElements
                        displ_dir(:,geo.Elements(k,geo.MasterFaces(4,:))) = true;
                    end
                case 3
                    % Fix broad end of TA
                    displ_dir(:,geo.Elements(8,geo.MasterFaces(4,:))) = true;
                    displ_dir(:,geo.Elements(8,geo.MasterFaces(2,:))) = true;
                    % For quick release test we also constrain the movement
                    % of the loose end to the y axis
                    if ~this.ICCompMode
                        displ_dir([1 3],geo.Elements(6,geo.MasterFaces(3,:))) = true;
                    end
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
                        len = box(2)-box(1);
                        velo_dir(1,geo.Elements(2,geo.MasterFaces(2,:))) = true;
                        totaldistance = (f.lambdafopt-1)*len;
                        velo_dir_val(velo_dir) = totaldistance/this.icMovetime;
                    case 2
                        for k = 1:4
                           velo_dir(2,geo.Elements(k,geo.MasterFaces(3,:))) = true;
                        end
                        totaldistance = (f.lambdafopt-1)*(box(4)-box(3));
                        velo_dir_val(velo_dir) = -totaldistance/this.icMovetime;
                    case 3
                        len = box(4)-box(3);
                        % Set only for y direction to nonzero value
                        velo_dir(2,geo.Elements(6,geo.MasterFaces(3,:))) = true;
                        totaldistance = (f.lambdafopt-1)*len;
                        velo_dir_val(velo_dir) = -totaldistance/this.icMovetime;
                        % Set zero for rest
                        velo_dir(:,geo.Elements(6,geo.MasterFaces(3,:))) = true;
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
            x0file = fullfile(QuickRelease.OutputDir,sprintf('geo%d_ic.mat',geonr));
            if exist(x0file,'file') ~= 2
                mc = QuickRelease(geonr, true);
                m = muscle.Model(mc);
                mu = m.DefaultMu;
                mu(1) = .001;
                [t, y] = m.simulate(mu);
%                 m.plot(t,y);
                yfull = m.System.includeDirichletValues(t,y);
                x0 = yfull(:,end);%#ok
                save(x0file,'x0','t','y','m');
                fprintf('Initial conditions file created. Please re-run.\n');
                return;
            end
            
            file = fullfile(QuickRelease.OutputDir,sprintf('geo%d.mat',geonr));
            saveit = true;
            if exist(file,'file') == 2
                load(file);
                saveit = false;
            else
                mc = QuickRelease(geonr, false);
                m = muscle.Model(mc);
            end
%             geo = mc.PosFE.Geometry;
            
            pm = PlotManager(false,2,2);
%             pm = PlotManager;
            pm.AutoTickMarks = 0;
            pm.LeaveOpen = true;
            
            c = ColorMapCreator;
            c.useJet([0.01 0.05 0.1 .5]);
            
            mus = repmat(m.DefaultMu,1,2);
            mus(1,:) = [.1 1];
            m.Data.ParamSamples = mus;
            nparams = length(mus);
            for k=1:nparams
                mu = mus(:,k);
                for inidx = 1:m.System.InputCount
                    [t, y] = m.simulate(mu,inidx);
                    o = m.Config.getOutputOfInterest(m, t, y);
                    
                    % Remove the first two entries as the activation is
                    % increased there
%                     o = o(:,3:end);

                    if saveit
                        save(file,'m','y','o','t','mc');
                    end
                    
                    h = pm.nextPlot(sprintf('force_velo_load%g_visc%g',mc.Loads(inidx),mu(1)),...
                        sprintf('Force / velocity plot\nLoad: %g[g] (eff. pressure %g[kPa]), viscosity:%g [mNs/m]',mc.Loads(inidx),mc.Pressures(inidx),mu(1)),...
                        'Velocity [mm/ms]','Force [N]');
                    plot(h,o(1,:),o(2,:),'r');
                    
%                     h = pm.nextPlot(sprintf('velo_visc%g',v),...
%                         sprintf('Averaged X-velocity [mm/ms] of right face for different rates and activation level\nViscosity=%g',v),...
%                         'Activation increase rate [1/ms]','alpha [-]');
%                     surf(h,1./ramptimes, alphavals, abs(pos(:,:,2))','FaceColor','interp','EdgeColor',ec);
%                     set(h,'XScale','log');
%                     view(-136, 56);
                end
            end
            
%             m.plotGeometrySetup(pm);
%             
            pm.done;
%             pm.FilePrefix = sprintf('case1_geo%d',geonr);
%             pm.ExportDPI = 200;
%             pm.SaveFormats = {'jpg'};
%             pm.savePlots(QuasiStaticTest.OutputDir,'Close',true);
        end
       
    end
    
end

