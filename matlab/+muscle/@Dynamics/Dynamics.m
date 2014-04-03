classdef Dynamics < dscomponents.ACompEvalCoreFun
    
    
    properties
       System;
       
       c10 = 6.352e-10;
       c01 = 3.627;
    end
    
    methods
        function this = Dynamics(sys)
            this.System = sys;
            this.MultiArgumentEvaluations = false;
            
            dfe = sys.DisplFE;
            
            dirvals = length(sys.bc_dir_val);
            d = dfe.NumNodes * 6 - dirvals;
            this.xDim = d;
            this.fDim = d;
        end
        
        function duv = evaluateCoreFun(this, uv, t, mu)
            sys = this.System;
            g = sys.Model.Geometry;
            dfe = sys.DisplFE;
            globidx = sys.globidx;
            
            fielddofs = dfe.NumNodes*3;
            
            % Include dirichlet values to state vector
            uvall = zeros(2*fielddofs,1);
            uvall(sys.dof_idx) = uv;
            uvall(sys.bc_dir_idx) = sys.bc_dir_val;
            
            % Init dy to yall vector.
            duv = zeros(size(uvall));
            % THIS ALREADY FULFILLS THE u' = v ODE PART!
            duv(1:fielddofs) = uvall(fielddofs+1:end);
            
            if t > 0
                a  =5 ;
            end
            
            dpe = dfe.DofsPerElement;
            ng = g.NumGaussp;
            ne = dfe.NumElems;
            for m = 1:ne
                %elem = dfe.elems(m,:);
                elemidx = globidx(:,:,m);
                
                integrand = zeros(3,dpe);
                for gp = 1:ng
                    
                    % TODO evaluate pressure at gp
                    p = .1;
                    
                    pos = 3*(gp-1)+1:3*gp;
                    dtn = dfe.transgrad(:,pos,m);
                    % Get coefficients for nodes of current element
                    c = uvall(elemidx);
                    
                    F = c * dtn;
                    
                    % Invariant I1
                    I1 = sum(sum((c'*c) .* (dtn*dtn')));
                    
                    P = p*inv(F)' + 2*(this.c10 + I1*this.c01)*F;
                    
                    integrand = integrand + g.gaussw(gp) * P * dtn' * dfe.elem_detjac(m,gp);
                    
                end
                % We have v' + K(u) = 0, so the values of K(u) must be
                % written at the according locations of v'; those are the
                % same but +fielddofs later each
                outidx = elemidx + fielddofs;
                duv(outidx) = duv(outidx) + integrand;
            end
            % Remove values at dirichlet nodes
            %dy(sys.bc_dir_idx)
            duv(sys.bc_dir_idx) = [];
        end
        
%         function dy = evaluateCoreFun(this, y, t, mu)%#ok
%             error('evaluate is overridden directly.');
%         end
        
%         function J = getStateJacobian(this, y, t, mu)
%             
%         end
        
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
        function fx = evaluateComponents(this, pts, ends, ~, ~, x, ~, mu)
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
        
        function dfx = evaluateComponentPartialDerivatives(this, pts, ends, idx, deriv, self, x, t, mu, dfxsel)%#ok
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
