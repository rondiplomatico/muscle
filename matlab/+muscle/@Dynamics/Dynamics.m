classdef Dynamics < dscomponents.ACompEvalCoreFun
    
    
    properties
       
       c10 = 6.352e-10;
       c01 = 3.627;
    end
    
    methods
        function this = Dynamics(sys)
            this = this@dscomponents.ACompEvalCoreFun(sys);
            
            dfe = sys.DisplFE;
            
            dirvals = length(sys.bc_dir_val);
            d = dfe.NumNodes * 6 - dirvals + sys.PressureFE.NumNodes;
            this.xDim = d;
            this.fDim = d;
        end
        
        function evaluateCoreFun(varargin)
            error('Custom projection is implemented and evaluate overridden directly');
        end
        
        function duvw = evaluate(this, uvwdof, t)
            if any(isnan(uvwdof(:)))
                keyboard;
            end
            sys = this.System;
            g = sys.Model.Geometry;
            fe_displ = sys.DisplFE;
            fe_press = sys.PressureFE;
            globidx_disp = sys.globidx_displ;
            globidx_press = sys.globidx_pressure;
            
            dofs_displ = fe_displ.NumNodes*3;
            
            % Include dirichlet values to state vector
            uvwcomplete = zeros(2*dofs_displ + fe_press.NumNodes,1);
            uvwcomplete(sys.dof_idx) = uvwdof;
            uvwcomplete(sys.bc_dir_idx) = sys.bc_dir_val;
            
            % Init duv
            duvw = zeros(size(uvwcomplete));
            % THIS ALREADY FULFILLS THE u' = v ODE PART!
            duvw(1:dofs_displ) = uvwcomplete(dofs_displ+1:2*dofs_displ);
            
            dofsperelem_displ = fe_displ.DofsPerElement;
            dofsperelem_press = fe_press.DofsPerElement;
            num_gausspoints = g.NumGaussp;
            num_elements = fe_displ.NumElems;
            for m = 1:num_elements
                elemidx_displ = globidx_disp(:,:,m);
                elemidx_pressure = globidx_press(:,m);
                
                integrand_displ = zeros(3,dofsperelem_displ);
                integrand_press = zeros(dofsperelem_press,1);
                for gp = 1:num_gausspoints
                    
                    % Evaluate the pressure at gauss points
                    w = uvwcomplete(elemidx_pressure);
                    p = w' * fe_press.Ngp(:,gp,m);
                    
                    pos = 3*(gp-1)+1:3*gp;
                    dtn = fe_displ.transgrad(:,pos,m);
                    
                    % Get coefficients for nodes of current element
                    u = uvwcomplete(elemidx_displ);
                    % Deformation gradient
                    F = u * dtn;
                    
                    % Invariant I1
                    I1 = sum(sum((u'*u) .* (dtn*dtn')));
                    
%                     if any(any(isnan(inv(F))))
%                         keyboard;
%                     end
                    
                    P = p*inv(F)' + 2*(this.c10 + I1*this.c01)*F - 2*this.c01*(F*F')*F;
%                     if norm(P) > 1e-5
%                         keyboard;
%                     end
                    
                    weight = g.gaussw(gp) * fe_displ.elem_detjac(m,gp);
                    
                    integrand_displ = integrand_displ + weight * P * dtn';
                    
                    integrand_press = integrand_press + weight * (det(F)-1) * fe_press.Ngp(:,gp,m);
                end
                % We have v' + K(u) = 0, so the values of K(u) must be
                % written at the according locations of v'; those are the
                % same but +fielddofs later each
                outidx = elemidx_displ + dofs_displ;
                duvw(outidx) = duvw(outidx) + integrand_displ;
                
                % Update pressure value at according positions
                duvw(elemidx_pressure) = duvw(elemidx_pressure) + integrand_press;
            end
            % Remove values at dirichlet nodes
            duvw(sys.bc_dir_idx) = [];
            
            if any(isnan(duvw(:)))
                keyboard;
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
