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
        function this = Model(conf)
            if nargin < 1
                conf = muscle.DebugConfig;
            end
            % Creates a new muscle model
            this = this@models.BaseFullModel;
            this.Name = sprintf('FEM Muscle model %s',class(conf));
            
            this.SaveTag = sprintf('musclemodel_%s',class(conf));
            this.Data = data.ModelData(this);
            %this.Data.useFileTrajectoryData;
            
            this.System = muscle.System(this);
            
            % This defines a default behaviour for all muscle models.
            % Override in AModelConfig.configureModel for dependent
            % behaviour.
            this.DefaultMu = [1; 50; 0; 0];
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
%             this.System.f.test_Jacobian;
%             chk = this.Config.PosFE.test_JacobiansDefaultGeo;
% %             chk = chk && this.Config.PosFE.test_QuadraticBasisFun;
%             chk = chk && this.Config.PressFE.test_JacobiansDefaultGeo;
%             if ~chk
%                 error('Health tests failed!');
%             end
        end
        
        function [t, x, time, cache] = computeTrajectory(this, mu, inputidx)
            % Allows to also call prepareSimulation for any quantities set
            % by the AModelConfig class.
            this.Config.prepareSimulation(mu, inputidx);
            [t, x, time, cache] = computeTrajectory@models.BaseFullModel(this, mu, inputidx);
        end
        
        function t = getConfigTable(this)
            f = this.System.f;
            t = PrintTable('Configuration of Model %s',this.Name);
            t.HasHeader = true;
            t.addRow('\rho_0', 'c_{10} [kPa]','c_{01} [kPa]','b_1 [kPa]','d_1 [-]','P_{max} [kPa]','\lambda_f^{opt}');
            t.addRow(this.MuscleDensity, f.c10,f.c01,f.b1,f.d1,f.Pmax,f.lambdafopt);
            t.Format = 'tex';
        end
        
        function plotForceLengthCurve(this, pm)
            if nargin < 2
                pm = PlotManager(false,1,2);
                pm.LeaveOpen = true;
            end
            f = this.System.f;
            
            lambda = .5:.005:2;
            fl = (f.Pmax / f.lambdafopt) * f.ForceLengthFun(lambda/f.lambdafopt);
            markertf = max(0,(f.b1./lambda.^2).*(lambda.^f.d1-1));
%             markertf = (f.b1./lambda.^2).*(lambda.^f.d1-1);
            
            pos = find(markertf > max(fl)*1.4,1,'first');
            if ~isempty(pos)
                lambda = lambda(1:pos);
                fl = fl(1:pos);
                markertf = markertf(1:pos);
            end
            
            h = pm.nextPlot('force_length',sprintf('Force-Length curve for model %s',this.Name),'\lambda [-]','pressure [kPa]');
            plot(h,lambda,fl,'r',lambda,markertf,'g',lambda,fl+markertf,'b');
%             axis(h,[0 2 0 2]);
            legend(h,'Active','Passive','Total','Location','NorthWest');
            
            dfl = (f.Pmax / f.lambdafopt) * f.ForceLengthFunDeriv(lambda/f.lambdafopt);
            dmarkertf = (lambda>=1).*(f.b1./lambda.^3).*((f.d1-2)*lambda.^f.d1 + 2);
            h = pm.nextPlot('force_length_deriv',sprintf('Derivative of Force-Length curve for model %s',this.Name),'\lambda','deriv [kPa/ms]');
            plot(h,lambda,dfl,'r',lambda,dmarkertf,'g',lambda,dfl + dmarkertf,'b');
%             axis(h,[0 2 -7 9]);
            
            if this.System.UseCrossFibreStiffness
                markertf = max(0,(f.b1cf./lambda.^2).*(lambda.^f.d1cf-1));
                h = pm.nextPlot('force_length_xfibre',sprintf('Force-Length curve in cross-fibre direction for model %s',this.Name),'\lambda [-]','pressure [kPa]');
                plot(h,lambda,markertf,'r');
                axis(h,[0 2 0 150]);
                legend(h,'Passive cross-fibre pressure','Location','NorthWest');
                
                dmarkertf = (lambda >= 1).*(f.b1cf./lambda.^3).*((f.d1cf-2)*lambda.^f.d1cf + 2);
                h = pm.nextPlot('force_length_xfibre_deriv',sprintf('Derivative of Force-Length curve in cross-fibre direction for model %s',this.Name),'\lambda [-]','pressure [kPa]');
                plot(h,lambda,dmarkertf,'r');
            end
            if nargin < 2
                pm.done;
            end
        end
        
        function plotAnisotropicPressure(this)
%             pm = PlotManager(false,1,2);
            pm = PlotManager;
            pm.LeaveOpen = true;
            f = this.System.f;
            
            [lambda, alpha] = meshgrid(.02:.02:1.5,0:.01:1);
            
            fl = f.ForceLengthFun(lambda/f.lambdafopt);
            active = f.Pmax./lambda.*fl.*alpha;% + 1*max(0,(lambda-1)).*alpha;
            passive = max(0,(f.b1./lambda.^2).*(lambda.^f.d1-1));
            
            h = pm.nextPlot('aniso_pressure',sprintf('Pressure in fibre direction for model %s',this.Name),'Stretch \lambda','Activation \alpha');
            surf(h,lambda,alpha,active+passive,'EdgeColor','k','FaceColor','interp');
            
%             zlim(h,[0,f.Pmax*1.5]);
            
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
            if ~isempty(varargin) && isa(varargin{1},'PlotManager')
                varargin = [{'PM'} varargin];
            end
            x0 = this.System.x0.evaluate(1);
            [~, nf] = this.getResidualForces(0, x0);
            if ~isempty(nf)
                varargin(end+1:end+2) = {'NF',nf};
            end
            [varargout{1:nargout}] = this.System.plot(0,x0,varargin{:});
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
        
        function idx = getDirichletBCFaceIdx(this, elem, face, dim)
            if nargin < 4
                dim = 1:3;
            end
            geo = this.Config.PosFE.Geometry;
            idx_face = false(size(this.System.bool_u_bc_nodes));
            idx_face(dim,geo.Elements(elem,geo.MasterFaces(face,:))) = true;
            fidx = find(this.System.bool_u_bc_nodes & idx_face);
            [~, idx] = intersect(this.System.idx_u_bc_glob, fidx);
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