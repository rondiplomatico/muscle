classdef Dynamics < dscomponents.ACompEvalCoreFun
    
    properties
       c10 = 6.352e-10; % [kPa]
       c01 = 3.627; % [kPa]
       b1 = 2.756e-5; % [kPa]
       d1 = 43.373; % [-]
       Pmax = 7.3; % N/cm²
       lambdafopt = 1.2; % [-]
       alpha = .1;
    end
    
    methods
        function this = Dynamics(sys)
            this = this@dscomponents.ACompEvalCoreFun(sys);
            
            mc = sys.Model.Config;
            dfe = mc.PosFE;
            
            dirvals = length(sys.bc_dir_val);
            d = dfe.NumNodes * 6 - dirvals + mc.PressFE.NumNodes;
            this.xDim = d;
            this.fDim = d;
            this.computeSparsityPattern;
        end
        
        function evaluateCoreFun(varargin)
            error('Custom projection is implemented and evaluate overridden directly');
        end
        
        function computeSparsityPattern(this)
            sys = this.System;
            mc = sys.Model.Config;
            g = mc.Geometry;
            fe_pos = mc.PosFE;
            fe_press = mc.PressFE;
            
            N = fe_pos.NumNodes;
            M = fe_press.NumNodes;
            
            %% -I part in u'(t) = -v(t)
            i = (1:3*N)';
            j = ((1:3*N)+3*N)';
            
            globidx_disp = sys.globidx_displ;
            globidx_press = sys.globidx_pressure;
            
            dofs_displ = N*3;
            
            dofsperelem_displ = fe_pos.DofsPerElement;
            dofsperelem_press = fe_press.DofsPerElement;
            num_gausspoints = g.NumGaussp;
            num_elements = fe_pos.NumElems;
            for m = 1:num_elements
                elemidx_displ = globidx_disp(:,:,m);
                elemidx_velo = elemidx_displ + dofs_displ;
                elemidx_pressure = globidx_press(:,m);
                inew = elemidx_velo(:);
                one = ones(size(inew));
                for gp = 1:num_gausspoints
                    for k = 1:dofsperelem_displ
                        %% Grad_u K(u,v,w)
                        % xdim
                        i = [i; inew]; %#ok<*AGROW>
                        j = [j; one*elemidx_displ(1,k)];
                        
                        % ydim
                        i = [i; inew]; 
                        j = [j; one*elemidx_displ(2,k)]; 
                        
                        % zdim
                        i = [i; inew]; 
                        j = [j; one*elemidx_displ(3,k)]; 
                        
                        %% grad u g(u)
                        % dx
                        i = [i; elemidx_pressure(:)];
                        j = [j; ones(dofsperelem_press,1)*elemidx_displ(1,k)]; 
                        % dy
                        i = [i; elemidx_pressure(:)];
                        j = [j; ones(dofsperelem_press,1)*elemidx_displ(2,k)]; 
                        %dz
                        i = [i; elemidx_pressure(:)];
                        j = [j; ones(dofsperelem_press,1)*elemidx_displ(3,k)]; 
                    end
                    %% Grad_w K(u,v,w)
                    inew = elemidx_velo(:);
                    for k = 1:dofsperelem_press
                        i = [i; inew];
                        j = [j; ones(3*dofsperelem_displ,1)*elemidx_pressure(k)]; 
                    end
                end
            end
            J = sparse(i,j,ones(size(i)),6*N+M,6*N+M);
            % Remove values at dirichlet nodes
            J(:,sys.bc_dir_idx) = [];
            J(sys.bc_dir_idx,:) = [];
            
            this.JSparsityPattern = logical(J);
        end
        
        function res = test_Jacobian(this, varargin)
            % Overrides the random argument jacobian test as restrictions
            % on the possible x values (detF = 1) hold.
            
            %res = test_Jacobian@dscomponents.ACoreFun(this, varargin{:});
            x0 = this.System.x0.evaluate([]);
            res = test_Jacobian@dscomponents.ACoreFun(this, x0);
            
            % Check if sparsity pattern and jacobian matrices match
            Jp = this.JSparsityPattern;
            Jeff = Jp;
            Jeff(:) = false;
            J = this.getStateJacobian(x0);
            Jeff(logical(J)) = true;
            check = (Jp | Jeff) & ~Jp;
            if any(check(:))
                error('boo');
            end
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
    
%     methods(Access=private)
%         function gv = g(this, lambdafsq)
%             alpha = .1;
%             lambdaf = sqrt(lambdafsq);
%             gv = (this.b1/lambdafsq)*(lambdaf^this.d1-1);
%             ratio = lambdaf/this.lambdafopt;
%             fl = (-6.25*ratio*ratio + 12.5*ratio - 5.25) * (ratio >= .6) * (ratio < 1.4);
%             gv = gv + (this.Pmax/lambdaf)*fl*alpha;
%         end
%     end
 
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
