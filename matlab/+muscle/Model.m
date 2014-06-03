classdef Model < models.BaseFullModel
    % Model: Model for a FEM-discretized muscle model
    %
    % The global time unit for this model is milliseconds [ms] and the
    % spatial quantities are in [mm].
    % This results in pressure values of [kPa] and the forces K(u,v,w) are
    % measured in [mN] (milliNewton).
    %
    % @author Daniel Wirtz @date 2012-11-22
    
    properties(SetAccess=private)
        % Seed that can be used by random number generator instances in order to enable result
        % reproduction.
        % @type int @default 1
        RandSeed = 1;
        
        Config;
        
        MuscleDensity = 1.1e-6; % [kg/mm³] (1100kg/m³)
    end
    
    methods
        function this = Model(conf)
            if nargin < 1
                conf = muscle.DebugConfig;
            end
            % Creates a new muscle model
            this = this@models.BaseFullModel;
            this.Name = sprintf('FEM Muscle model%s','');
            
            this.SaveTag = sprintf('musclemodel_%s','');
            this.Data = data.ModelData(this);
            %this.Data.useFileTrajectoryData;
            
            this.System = muscle.System(this);
            
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
        
        function plotForceLengthCurve(this)
            pm = PlotManager(false,1,2);
            pm.LeaveOpen = true;
            f = this.System.f;
            h = pm.nextPlot('force_length',sprintf('Force-Length curve for model %s',this.Name),'\lambda/\lambda_{opt}','force');
            ratio = 0:.01:2;
            plot(h,ratio,f.ForceLengthFun(ratio),'r');
            h = pm.nextPlot('force_length_deriv',sprintf('Derivative of Force-Length curve for model %s',this.Name),'\lambda/\lambda_{opt}','deriv');
            plot(h,ratio,f.ForceLengthFunDeriv(ratio),'b');
            pm.done;
        end
        
        function plotGeometrySetup(this, pm)
            args = {};
            if nargin == 2
                args = {'PM',pm};
            end
            x0 = this.System.x0.evaluate(1);
            [~, nf] = this.getResidualForces(0, x0);
            if ~isempty(nf)
                args(end+1:end+2) = {'NF',nf};
            end
            this.System.plot(0,x0,args{:});
        end
        
        function [residuals_dirichlet, residuals_neumann] = getResidualForces(this, t, uvw)
            sys = this.System;
            posvals = length(sys.bc_dir_displ_idx);
            residuals_dirichlet = zeros(posvals+length(sys.bc_dir_velo_idx),length(t));
            residuals_neumann = zeros(length(sys.bc_neum_forces_nodeidx),length(t));
            pos_dofs = this.Config.PosFE.Geometry.NumNodes * 3;
            dyall = zeros(2*pos_dofs + this.Config.PressFE.Geometry.NumNodes,1);
            for k=1:length(t)
                dy = sys.f.evaluate(uvw(:,k),t(k));
                % Place in global vector
                dyall(sys.dof_idx_global,:) = dy;
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
            idx_face = false(size(this.System.bc_dir_displ));
            idx_face(dim,geo.Elements(elem,geo.MasterFaces(face,:))) = true;
            fidx = find(this.System.bc_dir_displ & idx_face);
            [~, idx] = intersect(this.System.bc_dir_displ_idx, fidx);
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
    end
    
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
                this = loadobj@models.BaseFullModel(this, sobj);
            else
                this = loadobj@models.BaseFullModel(this);
            end
        end
    end
end