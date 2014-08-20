classdef Dynamics < dscomponents.ACompEvalCoreFun
    % This class implements the nonlinear continuum mechanics as described
    % in @cite Heidlauf2013 .
    
    properties
       %% Isotropic parameters (Moonley-Rivlin)
       c10 = 6.352e-10; % [kPa]
       c01 = 3.627; % [kPa]
       
       %% Anisotropic parameters (Markert)
       b1 = 2.756e-5; % [kPa]
       d1 = 43.373; % [-]
       
       % Cross-fibre markert part
       b1cf = 53163.72204148964; % [kPa] = [N/mm²]
       d1cf = 0.014991843974911; % [-]
       
       Pmax = 73; % [kPa], in Paper 7.3N/cm², but kPa = 0.1N/cm² 
       
       lambdafopt = 1.2; % [-]
       
       % The activation of the muscle at time t
       %
       % @type function_handle @default @(t)0
       alpha = @(t)0; % [-]
       
       % The force-length function as function handle
       %
       % This function describes the force-length relation for active
       % force.
       %
       % The alternative function using gaussian-type shapes is from
       % @cite Guenther2007
       %
       % @type function_handle @default Quadratic like in @cite
       % Heidlauf2013
       ForceLengthFun = @(ratio)(-6.25*ratio.*ratio + 12.5*ratio - 5.25) .* (ratio >= .6) .* (ratio <= 1.4);
       % Alternative
       % ForceLengthFun = @(ratio)(ratio<=1).*exp(-((1-ratio)/.57).^4) + (ratio>1).*exp(-((ratio-1)/.14).^3);
       
       % The derivative of the force-length function as function handle
       %
       % This function describes the derivative of the force-length
       % relation for active force.
       %
       % The alternative function using gaussian-type shapes is from
       % @cite Guenther2007
       %
       % @type function_handle @default Quadratic like in @cite
       % Heidlauf2013
       ForceLengthFunDeriv = @(ratio)(12.5*ratio.*(1-ratio)) .* (ratio >= .6) .* (ratio <= 1.4);
       % Alternative
       % ForceLengthFunDeriv = @(ratio)(ratio<=1).*((1/.57)*(((1-ratio)/.57).^3).*exp(-((1-ratio)/.57).^4)) ...
       % - (ratio > 1) .* ((1/.14) .* (((ratio-1)/.14).^2) .* exp(-((ratio-1)/.14).^3));
    end
    
    properties(SetAccess=private)
        APExp;
    end
    
    properties(Transient, SetAccess=private)
        % Prepared arguments for APExpansion
%         muprep;
        
        LastBCResiduals;
    end
    
    properties(Transient)
        % Helper value for QuickReleaseTests (or others) that use a
        % function handle with certain alpha ramp time for different
        % simulations. Used in getSimCacheExtra to uniquely identify a
        % simulation in the cache
        RampTime;
    end
    
    properties(Transient, Access=private)
        % Cached quantity from this.System.UseDirectMassInversion for
        % faster evaluation of dynamics.
        usemassinv;
        
        % Cached value for cross fibre computations (speed)
        crossfibres = false;
    end
    
    methods
        function this = Dynamics(sys)
            this = this@dscomponents.ACompEvalCoreFun(sys);
            
            %% Load AP expansion
%             d = fileparts(which('muscle.Dynamics'));
%             s = load(fullfile(d,'AP'));
%             s.kexp.Ma = s.kexp.Ma(1,:);
%             this.APExp = s.kexp;
        end
        
        function configUpdated(this)
            mc = this.System.Model.Config;
            if ~isempty(mc)
                geo = mc.PosFE.Geometry;

                dirvals = length(this.System.bc_dir_val);
                d = geo.NumNodes * 6 - dirvals + mc.PressFE.Geometry.NumNodes;
                this.xDim = d;
                this.fDim = d;
                this.computeSparsityPattern;
            end
        end
        
        function prepareSimulation(this, mu)
            prepareSimulation@dscomponents.ACompEvalCoreFun(this, mu);
            mc = this.System.Model.Config;
%             if ~isempty(mc.FibreTypes)
%                 this.muprep = [mc.FibreTypes; ones(size(mc.FibreTypes))*mu(2)];
%             end
            if ~isempty(mc.Pool)
                mc.Pool.prepareSimulation(this.System.Model.T,mu(2));
            end
            this.usemassinv = this.System.UseDirectMassInversion;
            this.crossfibres = this.System.UseCrossFibreStiffness;
        end
        
        function evaluateCoreFun(varargin)
            error('Custom projection is implemented and evaluate overridden directly');
        end
        
        function fx = evaluateComponentSet(this, nr, x, t)
            % Computes the full or reduced component functions of the given point set.
            %
            % Parameters:
            % nr: The number of the PointSet to use. @type integer
            % x: The state space location `\vx` @type colvec<double>
            % t: The corresponding times `t` for the state `\vx` @type double
            %
            % See also: PointSet
            if ~isempty(this.V)
                fx = this.evaluate(this.V*x,t);
            else
                fx = this.evaluate(x,t);
            end
            fx = fx(this.PointSets{nr},:);
        end

        function fx = evaluateComponentSetMulti(this, nr, x, t, mu)
            % Computes the full component functions of the given point set.
            %
            % Parameters:
            % nr: The number of the PointSet to use. @type integer
            % x: The state space locations `\vx` @type matrix<double>
            % t: The corresponding times `t` for the state `\vx` @type
            % rowvec<double>
            % mu: The corresponding parameters `\mu`. Can be a single
            % parameter or as many as the size of x @type matrix<double>
            %
            % See also: PointSet
            fx = this.evaluateMulti(x,t,mu);
            fx = fx(this.PointSets{nr},:);
        end
        
        function res = test_Jacobian(this, full)
            % Overrides the random argument jacobian test as restrictions
            % on the possible x values (detF = 1) hold.
            %
            % Currently the tests using viscosity are commented out as we
            % assume linear damping, which is extracted as extra `A(t,\mu)`
            % part in the models' system
            
            if nargin < 2
                full = false;
            end
            
            m = this.System.Model;
            if full
                [t,y] = m.simulate;
            else
                t = 100;
                y = this.System.x0.evaluate(this.System.Model.DefaultMu);
            end
            
            % Use nonzero t to have an effect
            res = test_Jacobian@dscomponents.ACoreFun(this, y, t, this.System.Model.DefaultMu);
            
            % Check if sparsity pattern and jacobian matrices match
%             Jp = this.JSparsityPattern;
%             Jeff = Jp;
%             for k=1:length(t)
%                 J = this.getStateJacobian(y(:,k),t(k));
%                 Jeff(:) = false;
%                 Jeff(logical(J)) = true;
%                 check = (Jp | Jeff) & ~Jp;
%                 res = res && ~any(check(:));
%             end
        end
        
        function copy = clone(this)
            % Create new instance
            copy = muscle.Dynamics(this.System);
            
            % Call superclass clone (for deep copy)
            copy = clone@dscomponents.ACompEvalCoreFun(this, copy);
            
            % Copy local properties
            % No local properties are to be copied here, as so far everything is done in the
            % constructor.
        end
    end
 
    methods(Access=protected)
        function fx = evaluateComponents(this, pts, ends, ~, ~, x, ~)
            % This is the template method that actually evaluates the components at given points
            % and values.
            %
            % @attention This method must be able to handle vector-arguments
            % for `\vx,t,\vmu`!
            %
            % Parameters:
            % pts: The components of `\vf` for which derivatives are required @type rowvec<integer>
            % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
            % indicates an input value relevant for the `i`-th point evaluation, i.e.
            % `f_i(\vx) = f_i(\vx(ends(i-1)+1{:}ends(i)));` @type rowvec<integer>
            % idx: The indices of `\vx`-entries in the global `\vx` vector w.r.t the `i`-th
            % point, e.g. `xglobal(i-1:i+1) = \vx(ends(i-1):ends(i))` @type rowvec<integer>
            % self: The positions in the `\vx` vector that correspond to the `i`-th output
            % dimension, if applicable (usually `f_i` depends on `x_i`, but not necessarily)
            % @type rowvec<integer>
            % x: A matrix `\vX` with the state space locations `\vx_i` in
            % its columns. In rows end(i-1)+1:end(i) it contains the states
            % relevant to evaluate the i-th component of �\vf�.
            % States occur multiply in �\vx� if different components
            % of �\vf� depend on these states.
            % @type matrix<double>
            % t: The corresponding times `t_i` for each state `\vx_i` @type rowvec<double>
            % mu: The corresponding parameters `\mu_i` for each state `\vx_i`, as column matrix
            % @type matrix<double>
            %
            % Return values:
            % fx: A matrix with pts-many component function evaluations `f_i(\vx)` as rows and as
            % many columns as `\vX` had.
            
           
        end
        
        function dfx = evaluateComponentPartialDerivatives(this, pts, ends, idx, deriv, self, x, t, dfxsel)%#ok
            % Computes specified partial derivatives of `f` of the components given by pts and
            % the selected partial derivatives by dfxsel.
            %
            % Parameters:
            % pts: The components of `f` for which derivatives are required @type
            % rowvec<integer>
            % ends: At the `i`-th entry it contains the last position in the `\vx` vector that
            % indicates an input value relevant for the `i`-th point evaluation, i.e.
            % `f_i(\vx) = f_i(\vx(ends(i-1){:}ends(i)));` @type rowvec<integer>
            % idx: The indices of `\vx`-entries in the global `\vx` vector w.r.t the `i`-th
            % point, e.g. `xglobal(i-1:i+1) = \vx(ends(i-1):ends(i))` @type rowvec<integer>
            % deriv: The indices within `\vx` that derivatives are required for.
            % @type rowvec<integer>
            % self: The positions in the `\vx` vector that correspond to the `i`-th output
            % dimension, if applicable (usually `f_i` depends on `x_i`, but not necessarily)
            % @type rowvec<integer>
            % x: The state space location `\vx` @type colvec<double>
            % t: The corresponding times `t` for the state `\vx` @type double
            % mu: The corresponding parameter `\mu` for the state `\vx` @type colvec<double>
            % dfxsel: A derivative selection matrix. Contains the mapping for each row of x to
            % the output points pts. As deriv might contain less than 'size(x,1)' values, use
            % 'dfxsel(:,deriv)' to select the mapping for the actually computed derivatives.
            %
            % Return values:
            % dfx: A column vector with 'numel(deriv)' rows containing the derivatives at all
            % specified pts i with respect to the coordinates given by 'idx(ends(i-1):ends(i))'
            %
            % See also: setPointSet
            
        end
    end    
end
