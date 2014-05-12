classdef Model < models.BaseFullModel
    % Model: Model for a FEM-discretized muscle model
    %
    % The global time unit for this model is milliseconds [ms] and the
    % spatial quantities are in [mm].
    % This results in pressure values of [kPa].
    %
    % @author Daniel Wirtz @date 2012-11-22
    
    properties
        % Seed that can be used by random number generator instances in order to enable result
        % reproduction.
        % @type int @default 1
        RandSeed = 1;
        
        Config;
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
            
            this.setConfig(conf);
            
            % True timestepping in ODE solver
            %   this.System.MaxTimestep = 1e-5; % [ms]
%             this.dt = 0.01;
%             this.T = 1;
%             this.ODESolver = solvers.ExplEuler;
            
            this.T = 10; % [ms]
            this.dt = .01; % [ms]
            s = solvers.MLode15i;
            % Use relatively coarse precision settings. This skips fine
            % oscillations but yields the correct long time results.
            s.RelTol = 1e-2;
            s.AbsTol = 1e-2;
            this.ODESolver = s;
            this.System.MaxTimestep = []; %model.dt;
            
            %% Health tests
%             this.System.f.test_Jacobian;
%             chk = this.Config.PosFE.test_JacobiansDefaultGeo;
% %             chk = chk && this.Config.PosFE.test_QuadraticBasisFun;
%             chk = chk && this.Config.PressFE.test_JacobiansDefaultGeo;
%             if ~chk
%                 error('Health tests failed!');
%             end
        end
        
        function setConfig(this, value)
            if ~isa(value, 'muscle.AModelConfig')
                error('Config must be a muscle.AModelConfig instance');
            end
            this.Config = value;
            this.System.configUpdated;
        end
        
        function varargout = plot(this, varargin)
            % plots some interesting states of the model
            %
            % See also: musclefibres.MuscleFibreSystem
            [varargout{1:nargout}] = this.System.plot(varargin{:});
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