classdef Model < models.BaseFullModel
    % Model: Model for a FEM-discretized muscle model
    %
    % The global time unit for this model is milliseconds [ms] and the
    % spatial quantities are in [mm].
    % This results in pressure values of [kPa] and the forces K(u,v,w) are
    % measured in [mN] (milliNewton).
    %
    % @author Daniel Wirtz @date 2012-11-22
    
    properties
        MuscleDensity = 1.1e-6; % [kg/mm³] (1100kg/m³)
        
        % Any arguments to add to any model.plot command for convenience.
        PlotterDefaultArgs = {};
    end
    
    properties(SetAccess=private)
        % Seed that can be used by random number generator instances in order to enable result
        % reproduction.
        % @type int @default 1
        RandSeed = 1;
        
        Config;
        
        Gravity = 9.80665; % [m/s²]
    end
    
    properties(Dependent)
        Geo;
    end
    
    methods
        function this = Model(conf, basedir)
            if nargin < 2
                basedir = KerMor.App.DataDirectory;
                if nargin < 1
                    conf = muscle.DebugConfig;
                end
            end
            % Creates a new muscle model
            name = sprintf('FEM Muscle model %s',class(conf));
            this = this@models.BaseFullModel(name);
            
            this.SaveTag = sprintf('musclemodel_%s',class(conf));
            this.Data = data.ModelData(this,basedir);
            %this.Data.useFileTrajectoryData;
            
            this.System = muscle.System(this);
            
            % This defines a default behaviour for all muscle models.
            % Override in AModelConfig.configureModel for dependent
            % behaviour.
            % Default is no activation ramp, no pressure and no input
            % current.
            this.DefaultMu = [1; 0; 0; 0
                % Anisotropic parameters muscle+tendon (Markert)
                2.756e-5; 43.373; 7.99; 16.6
                % Isotropic parameters muscle+tendon (Moonley-Rivlin)
                35.6; 3.86; 2310; 1.15e-3 % Micha, % 6.352e-10; 3.627 [Kpa] thomas alt 
                250; 1.2; .3]; 
            
            this.TrainingParams = [1 2];
            
            this.T = 10; % [ms]
            this.dt = .01; % [ms]
            s = solvers.MLode15i;
            % Use relatively coarse precision settings. This skips fine
            % oscillations but yields the correct long time results.
            s.RelTol = 1e-1;
            s.AbsTol = 1e-1;
            this.ODESolver = s;
            this.System.MaxTimestep = []; %model.dt;
            
            % Call the config-specific model configuration
            conf.Model = this;
            conf.configureModel(this);
            
            % Set the config to the model, triggering geometry related
            % pre-computations
            this.setConfig(conf);
            
            %% Health tests
            % Propagate the current default param
            this.System.prepareSimulation(this.DefaultMu, this.DefaultInput);
            
            fprintf('Running Jacobian health check..');
            res = this.System.f.test_Jacobian;
            fprintf('Done. Success=%d\n',res);
%             chk = this.Config.PosFE.test_JacobiansDefaultGeo;
% %             chk = chk && this.Config.PosFE.test_QuadraticBasisFun;
%             chk = chk && this.Config.PressFE.test_JacobiansDefaultGeo;
%             if ~chk
%                 error('Health tests failed!');
%             end
        end
        
        function [t,y] = simulateAndPlot(this, withResForce, varargin)
            if nargin < 2
                withResForce = true;
            end
            [t,y] = this.simulate(varargin{:});
            xargs = {};
            if (withResForce)
                [df,nf] = this.getResidualForces(t,y);
                xargs = {'NF',nf,'DF',df};
            end
            this.plot(t,y,xargs{:});
        end
        
        function [t, x, time, cache] = computeTrajectory(this, mu, inputidx)
            % Allows to also call prepareSimulation for any quantities set
            % by the AModelConfig class.
            this.Config.prepareSimulation(mu, inputidx);
            [t, x, time, cache] = computeTrajectory@models.BaseFullModel(this, mu, inputidx);
        end
        
        function t = getConfigTable(this, mu)
            if nargin < 2
                mu = this.DefaultMu;
            end
            f = this.System.f;
            t = PrintTable('Configuration of Model %s',this.Name);
            t.HasHeader = true;
            t.addRow('\rho_0', 'c_{10} [kPa]','c_{01} [kPa]','b_1 [kPa]','d_1 [-]','P_{max} [kPa]','\lambda_f^{opt}');
            t.addRow(this.MuscleDensity, mu(9),mu(10),mu(5),mu(6),mu(13),mu(14));
            t.Format = 'tex';
        end
        
        function plotForceLengthCurve(this, mu, pm)
            if nargin < 3
                pm = PlotManager(false,2,2);
                pm.LeaveOpen = true;
                if nargin < 2
                    mu = this.DefaultMu;
                end
            end
            f = this.System.f;
            f.setForceLengthFun(mu);
            
            markertfun = @(lam,b,d)max(0,(b./lam.^2).*(lam.^d-1));
            
            lambda = .2:.005:2;
            
            %% Plain Force-length function
            fl = f.ForceLengthFun(lambda/mu(14));
            h = pm.nextPlot('force_length_plain','Direct force-length curve of muscle/sarcomere',...
                sprintf('\\lambda/\\lambda_{opt} [-], \\lambda_{opt}=%g',mu(14)),...
                'force-length relation [-]');
            
            plot(h,lambda/mu(14),fl,'r');
            
            %% Effective force-length function for muscle tissue
            % effective signal from active part
            fl_eff = (mu(13)./lambda) .* fl;
            b1 = mu(5); d1 = mu(6);
            % Passive markert law
            markertf = markertfun(lambda,b1,d1);
            
            % Find a suitable position to stop plotting (otherwise the
            % passive part will steal the show)
            pos = find(markertf > max(fl_eff)*1.4,1,'first');
            if ~isempty(pos)
                lambda = lambda(1:pos);
                fl_eff = fl_eff(1:pos);
                markertf = markertf(1:pos);
            end
            
            h = pm.nextPlot('force_length_eff',...
                'Effective force-length curve of muscle material',...
                '\lambda [-]','pressure [kPa]');
            plot(h,lambda,fl_eff,'r',lambda,markertf,'g',lambda,fl_eff+markertf,'b');
            legend(h,'Active','Passive','Total','Location','NorthWest');
            
%             %% Effective force-length function derivative for muscle tissue
%             dfl = (mu(13)./lambda) .* f.ForceLengthFunDeriv(lambda/mu(14));
%             dmarkertf = (lambda>=1).*(b1./lambda.^3).*((d1-2)*lambda.^d1 + 2);
%             h = pm.nextPlot('force_length_eff_deriv',...
%                 'Effective force-length curve derivative of muscle material',...
%                 '\lambda','deriv [kPa/ms]');
%             plot(h,lambda,dfl,'r',lambda,dmarkertf,'g',lambda,dfl + dmarkertf,'b');
            
%             %% Cross fibre stuff
%             if this.System.UseCrossFibreStiffness
%                 error('fixme');
%                 markertf = max(0,(f.b1cf./lambda.^2).*(lambda.^f.d1cf-1));
%                 h = pm.nextPlot('force_length_xfibre',sprintf('Force-Length curve in cross-fibre direction for model %s',this.Name),'\lambda [-]','pressure [kPa]');
%                 plot(h,lambda,markertf,'r');
%                 axis(h,[0 2 0 150]);
%                 legend(h,'Passive cross-fibre pressure','Location','NorthWest');
%                 
%                 dmarkertf = (lambda >= 1).*(f.b1cf./lambda.^3).*((f.d1cf-2)*lambda.^f.d1cf + 2);
%                 h = pm.nextPlot('force_length_xfibre_deriv',sprintf('Derivative of Force-Length curve in cross-fibre direction for model %s',this.Name),'\lambda [-]','pressure [kPa]');
%                 plot(h,lambda,dmarkertf,'r');
%             end
            
            %% Passive force-length function for 100% tendon tissue
            % Passive markert law
            markertf = markertfun(lambda,mu(7),mu(8));
            h = pm.nextPlot('force_length_tendon',...
                'Effective force-length curve of tendon material (=passive)',...
                '\lambda [-]','pressure [kPa]');
            plot(h,lambda,markertf,'g');
            
            %% Effective force-length surface for muscle-tendon tissue
            pmap = this.System.MuscleTendonParamMapFun;
            % Sampled ratios
            tmr = 0:.03:1;
            %% Passive part
            % Resulting markert law params
            b1 = pmap(tmr,mu(5),mu(7));
            d1 = pmap(tmr,mu(6),mu(8));
            [LAM,TMR] = meshgrid(lambda,tmr);
            B = repmat(b1',1,length(lambda));
            D = repmat(d1',1,length(lambda));
            markertf = markertfun(LAM,B,D);
            
            %% Active part
            FL = (1-tmr)'*fl_eff;
            
            %% Plot dat stuff!
            h = pm.nextPlot('force_length_muscle_tendon',...
                'Effective force-length curve between muscle/tendon material',...
                '\lambda [-]','muscle/tendon ratio [m=0,t=1]');
            surfc(LAM,TMR,markertf+FL,'Parent',h,'EdgeColor','interp');
            zlabel('pressure [kPa]');
            zlim([0, 3*max(fl_eff(:))]);
                
            if nargin < 2
                pm.done;
            end
        end
        
        function plotAnisotropicPressure(this, mu)
            if nargin < 2
                mu = this.DefaultMu;
            end
%             pm = PlotManager(false,1,2);
            pm = PlotManager;
            pm.LeaveOpen = true;
            f = this.System.f;
            
            % Simply use muscle parameter values
            b1 = mu(1);
            d1 = mu(2);
            warning('Using muscle parameters only (ignoring tendon)');
            
            [lambda, alpha] = meshgrid(.02:.02:1.5,0:.01:1);
            
            fl = f.ForceLengthFun(lambda/mu(14));
            active = mu(13)./lambda.*fl.*alpha;% + 1*max(0,(lambda-1)).*alpha;
            passive = max(0,(b1./lambda.^2).*(lambda.^d1-1));
            
            h = pm.nextPlot('aniso_pressure',sprintf('Pressure in fibre direction for model %s',this.Name),'Stretch \lambda','Activation \alpha');
            surf(h,lambda,alpha,active+passive,'EdgeColor','k','FaceColor','interp');
            
            pm.done;
        end
        
        function plotActivation(this)
            pm = PlotManager;
            pm.LeaveOpen = true;
            f = this.System.f;
            
            h = pm.nextPlot('activation',sprintf('Activation curve for model %s',this.Name),'time [ms]','alpha [-]');
            plot(h,this.Times,f.alpha(this.scaledTimes),'r');
            pm.done;
        end
        
        function varargout = plotGeometrySetup(this, varargin)
            if isempty(this.System.mu)
                this.System.prepareSimulation(this.DefaultMu, this.DefaultInput);
            end
            if ~isempty(varargin) && isa(varargin{1},'PlotManager')
                varargin = [{'PM'} varargin];
            end
            x0 = this.System.x0.evaluate(this.System.mu);
            [~, nf] = this.getResidualForces(0, x0);
            if ~isempty(nf)
                varargin(end+1:end+2) = {'NF',nf};
            end
            [varargout{1:nargout}] = this.System.plot(0,x0,varargin{:});
        end
        
        function plotGeometryInfo(this, allnode, elemnr)
            if nargin < 3
                elemnr = 1;
                if nargin < 2
                    allnode = false;
                end
            end
            g = this.Config.PosFE.Geometry;
            g.plot(allnode,elemnr);
        end
        
        function [residuals_dirichlet, residuals_neumann] = getResidualForces(this, t, uvw)
            sys = this.System;
            num_bc = length(sys.idx_u_bc_glob)+length(sys.idx_v_bc_glob);
            residuals_dirichlet = zeros(num_bc,length(t));
            residuals_neumann = zeros(length(sys.bc_neum_forces_nodeidx),length(t));
            pos_dofs = this.Config.PosFE.Geometry.NumNodes * 3;
            dyall = zeros(2*pos_dofs + this.Config.PressFE.Geometry.NumNodes,1);
            for k=1:length(t)
                dy = sys.f.evaluate(uvw(:,k),t(k));
                % Place in global vector
                dyall(sys.idx_uv_dof_glob,:) = dy(1:sys.num_uvp_dof);
                % Select nodes that are exposed to neumann conditions (the
                % index is in global positions and not effective DoFs)
                residuals_neumann(:,k) = dyall(pos_dofs+sys.bc_neum_forces_nodeidx);
                
                residuals_dirichlet(:,k) = sys.f.LastBCResiduals;
            end
        end
        
        function idx = getFaceDofsGlobal(this, elem, faces, dim)
            % Returns the indices in the global uvw vector (including
            % dirichlet values) of the given faces in the given element.
            %
            % Parameters:
            % elem: The element index @type integer
            % faces: The faces whose indices to return. May also only be
            % one face. @type rowvec<integer>
            % dim: The dimensions which to select. @type rowvec<integer>
            % @default [1 2 3]
            %
            % Return values:
            % idx: The global indices of the face @type colvec<integer>
            if nargin < 4
                dim = 1:3;
            end
            geo = this.Config.PosFE.Geometry;
            idxXYZ = false(3,geo.NumNodes);
            for k=1:length(faces)
                idxXYZ(dim,geo.Elements(elem,geo.MasterFaces(faces(k),:))) = true;
            end
            idx = find(idxXYZ(:));
        end
        
        function idx = getPositionDirichletBCFaceIdx(this, elem, face, dim)
            % Returns the positions of dofs of a specified face 
            % within the boundary conditions residual vector.
            % Applies to dirichlet boundary conditions of POSITION (u)
            %
            % Parameters:
            % elem: The element the face belongs to @type integer
            % face: The face number of that element @type integer
            % dim: Optionally, specify a requested dimension to get the
            % indices for. Defaults to return all x,y,z components on every
            % node on the face. @type rowvec<integer> @default 1:3
            if nargin < 4
                dim = 1:3;
            end
            geo = this.Config.PosFE.Geometry;
            idx_face = false(size(this.System.bool_u_bc_nodes));
            idx_face(dim,geo.Elements(elem,geo.MasterFaces(face,:))) = true;
            fidx = find(this.System.bool_u_bc_nodes & idx_face);
            [~, idx] = intersect(this.System.idx_u_bc_glob, fidx);
        end
        
        function idx = getVelocityDirichletBCFaceIdx(this, elem, face, dim)
            % Returns the positions of dofs of a specified face 
            % within the boundary conditions residual vector.
            % Applies to dirichlet boundary conditions of VELOCITY (v)
            %
            % Parameters:
            % elem: The element the face belongs to @type integer
            % face: The face number of that element @type integer
            % dim: Optionally, specify a requested dimension to get the
            % indices for. Defaults to return all x,y,z components on every
            % node on the face. @type rowvec<integer> @default 1:3
            if nargin < 4
                dim = 1:3;
            end
            geo = this.Config.PosFE.Geometry;
            idx_face = false(size(this.System.bool_u_bc_nodes));
            idx_face(dim,geo.Elements(elem,geo.MasterFaces(face,:))) = true;
            fidx = find(this.System.bool_v_bc_nodes & idx_face);
            % Include the offset to the indices for velocity dofs
            fidx = fidx + geo.NumNodes*3;
            [~, idx] = intersect(this.System.idx_v_bc_glob, fidx);
            % Also include the offset of velocity components within the
            % dirichlet force vector
            idx = idx + length(this.System.val_u_bc);
        end
        
        function setConfig(this, value)
            if ~isa(value, 'muscle.AModelConfig')
                error('Config must be a muscle.AModelConfig instance');
            end
            this.Config = value;
            this.System.configUpdated;
        end
        
        function setGaussIntegrationRule(this, value)
            % Sets the gauss integration rule for the model.
            % 
            % See fem.BaseFEM for possible values. Currently 3,4,5 are
            % implemented.
            mc = this.Config;
            mc.PosFE.GaussPointRule = value;
            mc.PressFE.GaussPointRule = value;
            this.setConfig(mc);
        end
        
        function varargout = plot(this, varargin)
            % plots some interesting states of the model
            %
            % See also: musclefibres.MuscleFibreSystem
            [varargout{1:nargout}] = this.System.plot(varargin{:});
        end
        
        function value = get.Geo(this)
            value = this.Config.PosFE.Geometry;
        end
    end
    
%     methods(Access=protected)
%         function value = getSimCacheExtra(this)
%             % Return values:
%             % value: A column vector with additional values to distinguish
%             % the simulation from others (internal configurations) @type
%             % colvec<double>
%             value = this.System.f.RampTime;
%         end
%     end
    
    methods(Static)
        function test_ModelVersions
            % @TODO
            %
            % deformed reference geom
            % dirichlet nodes with not all 3 directions fixed, plotting
            
            
            m = muscle.Model(muscle.DebugConfig);
            mu = m.getRandomParam;
            [t,y] = m.simulate(mu);
            m.System.UseDirectMassInversion = true;
            [t,y] = m.simulate(mu);
            
            % Version with "constant" fibre activation forces
            m = muscle.Model(muscle.DebugConfig(2));
            [t,y] = m.simulate(mu);
            m.System.UseDirectMassInversion = true;
            [t,y] = m.simulate(mu);
            
            % Version with "neurophysiological" fibre activation forces
            m = muscle.Model(muscle.DebugConfig(3));
            [t,y] = m.simulate(mu);
            m.System.UseDirectMassInversion = true;
            [t,y] = m.simulate(mu);
            m.System.UseDirectMassInversion = false;
            f = m.System.f;
            f.Viscosity = 1;
            [t,y] = m.simulate(mu);
            m.System.UseDirectMassInversion = true;
            [t,y] = m.simulate(mu);
        end
        
        function test_JacobianApproxGaussRules(full)
            % Tests the precision of the analytical jacobian computation
            % using different gauss integration rules
            if nargin < 1
                full = false;
            end
            m = muscle.Model(muscle.DebugConfig(2));
            m.simulate(1);
            f = m.System.f;
            m.T = 1;
            m.dt = .2;
            f.test_Jacobian(full);
            m.setGaussIntegrationRule(4);
            f.test_Jacobian;
            m.setGaussIntegrationRule(5);
            f.test_Jacobian(full);
        end
    end
    
    methods(Static, Access=protected)
        function this = loadobj(this)
            if ~isa(this, 'muscle.Model')
                sobj = this;
                this = muscle.Model;
                this.RandSeed = sobj.RandSeed;
                this.Config = sobj.Config;
                this = loadobj@models.BaseFullModel(this, sobj);
            else
                this = loadobj@models.BaseFullModel(this);
            end
        end
    end
end