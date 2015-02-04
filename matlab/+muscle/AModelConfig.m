classdef AModelConfig < handle
    %AModelConfig
    
    properties
        Model;
    end
    
    properties(SetAccess=private)
        PosFE;
        
        PressFE;
    end
    
    properties(SetAccess=protected)
        FibreTypeWeights = [];
        
        Pool;
        
        % The coordinate system in which to interpret the applied pressure
        % of neumann boundary conditions.
        %
        % 'local' uses the normals on the faces as given in the reference
        % configuration, i.e. the "true" normals.
        %
        % 'global' uses the global coordinate system of the geometry, i.e.
        % the master element. This can be used to apply forces coming from
        % one fixed direction over a possibly noneven geometry surface.
        %
        % @type char @default 'local'
        NeumannCoordinateSystem = 'local';
        
        % The coordinate system in which to interpret the a0 vectors of
        % fibre directions.
        %
        % 'master' applies the given directions at the master element
        % coordinate system and transforms the directions according to the
        % transformation of the respective element.
        % 
        % 'reference' applies the given directions with respect to the
        % reference coordinate system, e.g. "as-is" at the gauss points.
        %
        % @type char @default 'master'
        a0CoordinateSystem = 'master';
    end
    
    properties
        % Determines the default value for maximum activation in activation
        % ramps.
        %
        % This is e.g. implicitly used by
        % muscle.Dynamics.prepareSimulation, when the alpha ramp is created
        % for positive mu(2) values (=ramp times).
        %
        % @type double @default 1
        ActivationRampMax = 1;
        
        % Determines the default number of milliseconds to wait before
        % activation is started.
        %
        % This is e.g. implicitly used by
        % muscle.Dynamics.prepareSimulation, when the alpha ramp is created
        % for positive mu(2) values (=ramp times).
        %
        % @type double @default 0
        ActivationRampOffset = 0;
        
        % Velocity conditions application function
        VelocityBCTimeFun;
    end
    
    methods
        function this = AModelConfig(geo)
            if isa(geo,'geometry.Cube8Node')
                pos_geo = geo.toCube27Node;
                press_geo = geo;
            elseif isa(geo,'geometry.Cube20Node') || isa(geo,'geometry.Cube27Node')
                pos_geo = geo;
                press_geo = geo.toCube8Node;
            else
                error('Scenario not yet implemented for geometry class "%s"', class(geo));
            end
            if isa(pos_geo,'geometry.Cube27Node')
                this.PosFE = fem.HexahedronTriquadratic(pos_geo);
            else
                this.PosFE = fem.HexahedronSerendipity(pos_geo);
            end
            this.PressFE = fem.HexahedronTrilinear(press_geo);
            %this.PressFE = fem.HexahedronSerendipity(press_geo.toCube20Node);
            %this.PressFE = fem.HexahedronTriquadratic(press_geo.toCube27Node);
        end
        
        function configureModel(this, model)
            % Overload this method to set model-specific quantities like
            % simulation time etc
            
            % Do nothing
        end
        
        function prepareSimulation(this, mu, inputidx)
            % Overload this method to initialize model-specific quantities
            % that are fixed for each simulation
            %
            % Called by override of computeTrajectory in muscle.Model
            
            % do nothing by default
        end
        
        function P = getBoundaryPressure(~, elemidx, faceidx)
            % Determines the neumann forces on the boundary.
            %
            % The unit for the applied quantities is kiloPascal [kPa]
            %
            % In the default implementation there are no force boundary
            % conditions.
            %
            % See also: NeumannCoordinateSystem
            P = [];
        end
        
        function u = getInputs(~)
            % Returns the inputs `u(t)` of the model.
            %
            % if neumann boundary conditions are used, this input is
            % multiplied with the mu(3) parameter, which determines the
            % maximum force that is applied. u(t) determines its temporal
            % strength.
            %
            % this.Model can be used to get access to the model this
            % configuration is applied to.
            %
            % Return values:
            % u: The cell array of input functions to use within this
            % model.
            %
            % @type cell @default {}
            u = {};
        end
        
        function x0 = getX0(~, x0)
            %% do nothing
        end
        
        function alpha = getAlphaRamp(this, ramptime, alphamax, starttime)
            % Creates a linearly increasing scalar function starting at
            % starttime milliseconds ranging from zero to alphamax over
            % ramptime.
            %
            % Parameters:
            % ramptime: The time over which to increase to alphamax. If
            % less or equal to zero, an all zero function is returned.
            % alphamax: The maximum value to achieve. @type double @default
            % AModelConfig.ActivationRampMax
            % starttime: The offset time (in milliseconds) to wait before
            % increasing the signal. @type double 
            % @default AModelConfig.ActivationRampOffset
            if ramptime <= 0
                alpha = @(t)0;
                return;
            end
            if nargin < 4
                starttime = this.ActivationRampOffset;
                if nargin < 3
                    alphamax = this.ActivationRampMax;
                end
            end
            alpha = @(t)(t >= starttime) .* (alphamax * (((t-starttime)<ramptime).*(t-starttime)/ramptime + (t>=ramptime+starttime)));
        end
        
        function tmr = getTendonMuscleRatio(~, ~)
            % Returns the [0,1] ratio between tendon and muscle at all
            % gauss points of all elements
            %
            % This method simply returns an empty ratio, meaning muscle only.
            %
            % Parameters:
            % x: A 3xn vector of coordinates at which to get the
            % tendonmuscle ratio @type matrix<double>
            %
            % Return values:
            % tmr: A row vector of tendonmuscle ratio values in [0 1] for
            % each of the n locations
            tmr = []; % zeros(1,size(x,2));
        end
    end
    
    methods(Access=protected)
        function anull = seta0(~, anull)
            % do nothing!
        end
        
        function [velo_dir, velo_dir_val] = setVelocityDirichletBC(~, velo_dir, velo_dir_val)
            % Determines the dirichlet velocities.
            %
            % The unit for the applied quantities is [mm/ms] = [m/s]
            %
            % In the default implementation there are no velocity
            % conditions.
        end        
        
        function ftw = getFibreTypeWeights(this)
            % This is a lazy pre-implementation as fullmuscle.Models
            % always have fibre types and thus weights.
            %
            % This method simply returns an all-zero weighting.
            fe = this.PosFE;
            geo = fe.Geometry;
            ftw = zeros(fe.GaussPointsPerElem,length(this.FibreTypes),geo.NumElements);
        end
    end
    
    methods(Abstract, Access=protected)
        displ_dir = setPositionDirichletBC(this, displ_dir);
    end
    
    methods(Sealed)
        function [displ_dir, velo_dir, velo_dir_val] = getBC(this)
            N = this.PosFE.Geometry.NumNodes;
            displ_dir = false(3,N);
            displ_dir = this.setPositionDirichletBC(displ_dir);
            velo_dir = false(3,N);
            velo_dir_val = zeros(3,N);
            [velo_dir, velo_dir_val] = this.setVelocityDirichletBC(velo_dir, velo_dir_val);
            
            if any(any(displ_dir & velo_dir)) 
                error('Cannot impose displacement and velocity dirichlet conditions on same DoF');
            end
        end
        
        function anull = geta0(this)
            fe = this.PosFE;
            g = fe.Geometry;
            anull = zeros(3,fe.GaussPointsPerElem,g.NumElements);
            anull = this.seta0(anull);
            % Normalize anull vectors
            for m = 1:g.NumElements
                anull(:,:,m) = anull(:,:,m) ./ ([1;1;1]*Norm.L2(anull(:,:,m)));
            end
        end
    end
    
end

