classdef SidePressureTest < experiments.AExperimentModelConfig
% Tests to investigate the difference between quasi-static and dynamic
% simulations.
%
% Tests can be run for three different geometries (Block, Belly & TA) and
% different viscosities
% 
%

    properties(Constant)
        % The time over which the muscle is activated
        ActivationTime = .5; % [ms]
        
        % The time waited at full activation before applying the load
        RelaxTime = 1; % [ms]
        
        % The time over which the load is increased/applied
        LoadRampTime = 4; % [ms]
    end

    properties(SetAccess=private)
        Loads;
        Pressures;
    end

    methods
        function this = SidePressureTest(varargin)
            this = this@experiments.AExperimentModelConfig(varargin{:});
            this.init;
            this.NeumannCoordinateSystem = 'global';
            this.ActivationRampMax = .5;
        end
        
        function configureModel(this, m)
            
            m.T = this.ActivationTime+this.RelaxTime+this.LoadRampTime+10;
            m.dt = m.T/200;
            m.DefaultInput = 1;
            m.DefaultMu(1) = .1;
            m.DefaultMu(2) = this.ActivationTime;
            % Apply boundary pressure (magnitude given by input!)
            m.DefaultMu(3) = 1;
            
            % Use the material set of
            f = m.System.f;
            
            m.DefaultMu(5) = 2.756e-5; % b1 [kPa]
            m.DefaultMu(6) = 43.373; % d1 [-]
            m.DefaultMu(9) = 6.352e-10; % c10 [kPa]
            m.DefaultMu(10) = 3.627; % c01 [kPa]
            
            % Cross-fibre markert part
            f.b1cf = 53163.72204148964/10; % [kPa] = [N/mm²]
            f.d1cf = 0.014991843974911; % [-]
            
            m.DefaultMu(13) = 250; % [kPa]
            m.DefaultMu(14) = 1; % [-]
            
            os = m.ODESolver;
            switch this.Options.GeoNr
                case 1
                    os.RelTol = .01;
                    os.AbsTol = .1;%.05;
                case {2,3}
                    os.RelTol = .01;
                    os.AbsTol = .5;
            end
            m.EnableTrajectoryCaching = true;
            m.System.UseCrossFibreStiffness = true;
        end
        
        function o = getOutputOfInterest(this, m, t, uvw)
%             geo = m.Geo;
            df = m.getResidualForces(t, uvw);
%             uvw = m.System.includeDirichletValues(t, uvw);
            switch this.Options.GeoNr
                case 1
                    idx = m.getPositionDirichletBCFaceIdx(1,1);
                    o(1,:) = mean(df(idx,:),1);
                case 2
                    xidx = m.getPositionDirichletBCFaceIdx(1:4,3,1);
                    yidx = m.getPositionDirichletBCFaceIdx(1:4,3,2);
                    zidx = m.getPositionDirichletBCFaceIdx(1:4,3,3);
                    vec = [mean(df(xidx,:),1); mean(df(yidx,:),1); mean(df(zidx,:),1)];
                    o(1,:) = Norm.L2(vec);
                    o(2:4,:) = vec;
                case 3
                    xidx = m.getPositionDirichletBCFaceIdx(6,3,1);
                    yidx = m.getPositionDirichletBCFaceIdx(6,3,2);
                    zidx = m.getPositionDirichletBCFaceIdx(6,3,3);
                    vec = [mean(df(xidx,:),1); mean(df(yidx,:),1); mean(df(zidx,:),1)];
                    o(1,:) = Norm.L2(vec);
                    o(2:4,:) = vec;
            end
        end
        
        function P = getBoundaryPressure(this, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is kiloPascal [kPa]
            %
            % See also: NeumannCoordinateSystem
            P = [];
            switch this.Options.GeoNr
                case 1
                    if elemidx == 2 && faceidx == 6
                        P = -1;
                    end
                case 2
                    if (elemidx == 9 || elemidx == 10) && faceidx == 6
                        P = -1;
                    end
                case 3
                    if (elemidx == 2 || elemidx == 3 ||...
                            elemidx == 9 || elemidx == 10) && faceidx == 6
                        P = -1;
                    end
            end
        end
        
        function u = getInputs(this)
            % loads in [g]
            loads = [0 100 200 500 1000];

            % convert to pressure:
            % [g]/1000 = [kg]
            % [kg]*[m/s²] = [N]
            % [N]*1000 = [mN]
            % [mN]/[mm²] = [kPa]
            switch this.Options.GeoNr
                case 1
                    a = this.PosFE.getFaceArea(2,6);
                case 2
                    a = this.PosFE.getFaceArea([11 12],[5 5]);
                case 3
                    a = this.PosFE.getFaceArea([2 3 9 10],[6 6 6 6]);
            end
            pressures = (loads/1000*this.Model.Gravity)*1000/a; % [kPa]
            u = {};
            start = this.ActivationTime+this.RelaxTime;
            % Configure a load input that increases over LoadRampTime
            % ms up to the set pressure
            for lidx = 1:length(loads)
                u{lidx} = this.getAlphaRamp(this.LoadRampTime,pressures(lidx),start);%#ok
            end
            this.Loads = loads;
            this.Pressures = pressures;
        end
    end
    
    methods(Access=protected)
        
        function geo = this.getGeometry(this)
            switch this.Options.GeoNr
                case 1
                    [pts, cubes] = geometry.Cube8Node.DemoGrid(0:20:60,[0 20],[0 20]);
                    geo = geometry.Cube8Node(pts, cubes);
                case 2
                    geo = Belly.getBelly(5, 50, 'Radius', 3.5, 'InnerRadius', 2.5, 'Gamma', 15);
                case 3
                    s = load(fullfile(fileparts(which(mfilename)),'..','CMISS','EntireTA.mat'));
                    geo = s.geo27;
            end
        end
        
        function displ_dir = setPositionDirichletBC(this, displ_dir)
            
            geo = this.PosFE.Geometry;
            %% Dirichlet conditions: Position (fix one side)
            % This is done for each test case.
            switch this.Options.GeoNr
                case 1
                    % Fix front
                    displ_dir(:,geo.Elements(1,geo.MasterFaces(1,:))) = true;
                    % Fix back
                    displ_dir(:,geo.Elements(3,geo.MasterFaces(2,:))) = true;
                    % Fix lower side in z direction
                    displ_dir(3,geo.Elements(2,geo.MasterFaces(5,:))) = true;
                case 2
                    % Both sides
                    for k = geo.NumElements-3:geo.NumElements
                        displ_dir(:,geo.Elements(k,geo.MasterFaces(4,:))) = true;
                    end
                    for k = 1:4
                        displ_dir(:,geo.Elements(k,geo.MasterFaces(3,:))) = true;
                    end
                    % and bottom
                    displ_dir(3,geo.Elements([11 12],geo.MasterFaces(5,:))) = true;
                case 3
                    % Fix broad end of TA
                    displ_dir(:,geo.Elements(8,geo.MasterFaces([2 4],:))) = true;
                    %displ_dir(:,geo.Elements(8,geo.MasterFaces(2,:))) = true;
                    % Fix thin end of TA
                    displ_dir(:,geo.Elements(6,geo.MasterFaces(3,:))) = true;
                    % Fix bottom of pressured area
                    displ_dir(3,geo.Elements([2 3 9 10],geo.MasterFaces(5,:))) = true;
            end
        end
        
        function anull = seta0(this, anull)
            switch this.Options.GeoNr
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
            
            file = fullfile(SidePressureTest.OutputDir,sprintf('geo%d.mat',geonr));
            saveit = true;
            if exist(file,'file') == 2
                load(file);
                saveit = false;
            else
                mc = SidePressureTest(geonr);
                m = muscle.Model(mc);    
            end
%             geo = mc.PosFE.Geometry;
            
            pm = PlotManager(false,2,2);
%             pm = PlotManager;
            pm.AutoTickMarks = 0;
            pm.LeaveOpen = true;
            
            c = ColorMapCreator;
            c.useJet([0.01 0.05 0.1 .5]);
            
            muv = [.01 .1 1];
            mus = repmat(m.DefaultMu,1,length(muv));
            mus(1,:) = muv;
            m.Data.ParamSamples = mus;
            nparams = size(mus,2);
            pi = ProcessIndicator('Running %d scenarios',nparams*m.System.InputCount,false,nparams*m.System.InputCount);
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
                        sprintf('Force plot\nLoad: %g[g] (eff. pressure %g[kPa]), viscosity:%g [mNs/m]',mc.Loads(inidx),mc.Pressures(inidx),mu(1)),...
                        'Time [ms]','Force [N]');
                    plot(h,t,o(1,:),'r','LineWidth',2);
                    
                    pi.step;
                end
            end
            pi.stop;
            
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

