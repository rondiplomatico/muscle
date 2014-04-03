classdef Model < models.BaseFullModel
    % Model: Model for a FEM-discretized muscle model
    %
    % The global time unit for this model is milliseconds [ms].
    %
    % @author Daniel Wirtz @date 2012-11-22
    
    properties
        % Seed that can be used by random number generator instances in order to enable result
        % reproduction.
        % @type int @default 1
        RandSeed = 1;
        
        Geometry;
    end
    
    methods
        function this = Model
            % Creates a new muscle model
            
            this.Name = sprintf('FEM Muscle model%s','');
            
            this.SaveTag = sprintf('musclemodel_%s','');
            this.Data = data.ModelData(this);
            %this.Data.useFileTrajectoryData;
            
            this.Geometry = cubegeom;
            this.System = muscle.System(this);
            
            % True timestepping in ODE solver
            %   this.System.MaxTimestep = 1e-5; % [ms]
            this.dt = 0.0001;
            this.T = 0.001;
            this.System.MaxTimestep = 0.0001; % [ms]
            this.ODESolver = solvers.ExplEuler;
            
%             this.T = 1; % [ms]
%             this.dt = .1; % [ms]
%             this.ODESolver = solvers.MLode15i;
%             this.System.MaxTimestep = []; %model.dt;

            
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