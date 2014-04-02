classdef fembase < handle
    %FEMBASE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        geo;
        
        N;
        
        gradN;
        
        % The mass matrix
        M;
        
        dofs;
        
        elems;
        
        elem_detjac;
        
        % transformed basis function gradients, stored in eldofs x gp*3
        % x element matrix (accessible in 20x3-chunks for each gauss point and element)
        transgrad;
    end
    
    properties(Dependent)
        NumDofs;
        NumElems;
        DofsPerElement;
    end
    
    methods
        function this = fembase(geo)
            if nargin < 1
                geo = cubegeom;
            end
            this.geo = geo;
        end
        
        function init(this)
            %% Compute edges from cubes.
            % Mind the corner numbering of the cubes! (see above)
            
            g = this.geo;
            el = this.elems;
            dof = this.dofs;
            np = size(dof,2);
            
%             %% Compute node to cube indices
%             pv = cell(np,1);
%             for i=1:np
%                 pv{i} = find(sum(el == i,2));
%             end
%             this.pts_cubes = pv;
            
            %% Precomputations
            % Iterate each cube
            gp = g.gaussp;
            % Number of gauss points
            ngp = size(gp,2);
            Mass = zeros(np,np);
            % Number of elements
            nel = size(el,1);
            % DOFs per Element
            eldofs = size(el,2);
            % Jacobian of element deformation at gauss points
            eljac = zeros(nel,ngp);
            % transformed basis function gradients, stored in eldofs x gp*3
            % x element matrix (accessible in 20x3-chunks for each gauss point and element)
            tg = zeros(eldofs,ngp*3,nel);
            %tg = zeros(eldofs,ngp*3,nel);
            % Iterate all volumes
            for m = 1:nel
                bval = zeros(eldofs,eldofs);
                elem = el(m,:);
                % Iterate all gauss points
                for gi = 1:ngp
                    % gauss point \xi
                    xi = gp(:,gi);
                    dNxi = this.gradN(xi);
                    % Jacobian of isogeometric mapping
                    Jac = this.dofs(:,elem)*dNxi;
                    %% Mass matrix related
                    % Jacobian at gauss point: get coordinates of dofs (8
                    % for trilinear, 20 for triquadratic) and multiply with
                    % gradient
                    eljac(m,gi) = det(Jac);
                    % Evaluate basis functions at xi
                    Nxi = this.N(xi);
                    % Add up [basis function values at xi] times [volume
                    % ratio at xi] times [gauss weight for xi]
                    bval = bval + g.gaussw(gi)*(Nxi*Nxi')*eljac(m,gi);
                    
                    %% Transformed Basis function gradients at xi
                    tg(:,3*(gi-1)+1:3*gi,m) = dNxi / Jac';
                end
                % Build up upper right part of symmetric mass matrix
                for j=1:eldofs
                    Mass(elem(j),elem(j:end)) = Mass(elem(j),elem(j:end)) + bval(j,j:end);
                end
            end
            this.transgrad = tg;
            this.M = sparse(Mass + Mass');
            this.elem_detjac = eljac;
        end
        
        function J = getJacobian(this, elem, xi)
            J = this.dofs(:,elem)*this.gradN(xi);
        end
        
        function plot(this, pm)
            if nargin < 2
                pm = PlotManager(false,2,2);
                pm.LeaveOpen = true;
            end
            this.geo.plot(pm);
            
            p = this.dofs;
            h = pm.nextPlot('dofs','Grid DOFs','x','y');
            plot3(h,p(1,:),p(2,:),p(3,:),'k.','MarkerSize',14);
            
            h = pm.nextPlot('mass','Mass matrix','node','node');
            mesh(h,full(this.M));
            
            pm.nextPlot('mass_pattern','Mass matrix pattern','dof','dof');
            spy(this.M);
            
            
            
            if nargin < 2
                pm.done;
            end
        end
    end
    
    methods
        function value = get.NumDofs(this)
            value = size(this.dofs,2);
        end
        
        function value = get.NumElems(this)
            value = size(this.elems,1);
        end
        
        function value = get.DofsPerElement(this)
            value = size(this.elems,2);
        end
    end
    
    methods(Static,Access=protected)
        function res = test_BasisFun(subclass)
            h = 1e-8;
            res = true;
            for k=1:10
                a = rand(3,1);
                A = repmat(a,1,3);
                diff = (subclass.N(A + h*eye(3)) - subclass.N(A))/h - subclass.gradN(a);
                res = res & all(diff(:) < 10*h);
            end
        end
    end
        
    methods(Static)
        function res = test_JacobiansDefaultGeo
            % Tests if the deformation jacobians using linear and quadratic
            % elements is the same for a default geometry
            ranges = {-1:1, -2:1, -2:2};
            res = true;
            for k = 1:length(ranges)
                [pts, cubes] = cubegeom.DemoCubeGrid(ranges{k},ranges{k});
                g = cubegeom(pts, cubes);
                tl = trilinear(g);
                tq = triquadratic(g);
                res = res && norm(tl.elem_detjac-tq.elem_detjac,'inf') < 1e-15;
            end
        end
    end
    
end

