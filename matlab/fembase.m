classdef fembase < handle
    %FEMBASE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Geometry;
        
        N;
        
        % Values of basis functions on gauss points
        Ngp;
        
        gradN;
        
        % The mass matrix
        M;
                
        elem_detjac;
        
        % transformed basis function gradients, stored in eldofs x gp*3
        % x element matrix (accessible in 20x3-chunks for each gauss point and element)
        transgrad;
    end
        
    methods
        function this = fembase(geometry)
            this.Geometry = geometry;
        end
        
        function init(this)
            %% Compute edges from cubes.
            % Mind the corner numbering of the cubes! (see above)
            
            g = this.Geometry;
            el = g.Elements;
            np = g.NumNodes;
            
            %% Precomputations
            % Iterate each cube
            gp = g.gaussp;
            % Number of gauss points
            ngp = g.GaussPointsPerElem;
            Mass = zeros(np,np);
            % Number of elements
            nel = size(el,1);
            % nodes per Element
            eldofs = size(el,2);
            
            % Basis function evaluation on Gauss points
            Ngpval = zeros(eldofs,ngp,nel);
            % Jacobian of element deformation at gauss points
            eljac = zeros(nel,ngp);
            % transformed basis function gradients, stored in eldofs x gp*3
            % x element matrix (accessible in 20x3-chunks for each gauss point and element)
            tg = zeros(eldofs,ngp*3,nel);
            % Iterate all volumes
            for m = 1:nel
                bval = zeros(eldofs,eldofs);
                elem = el(m,:);
                % Iterate all gauss points
                for gi = 1:ngp
                    % gauss point \xi
                    xi = gp(:,gi);
                    
                    %% N on gauss points
                    Ngpval(:,gi,m) = this.N(xi);
                    
                    %% Jacobian related
                    dNxi = this.gradN(xi);
                    % Jacobian of isogeometric mapping
                    Jac = g.Nodes(:,elem)*dNxi;
                    
                    %% Mass matrix related
                    % Jacobian at gauss point: get coordinates of nodes (8
                    % for trilinear, 20 for triquadratic) and multiply with
                    % gradient
                    eljac(m,gi) = det(Jac);
                    % Evaluate basis functions at xi
                    Nxi = this.N(xi);
                    % Add up [basis function values at xi] times [volume
                    % ratio at xi] times [gauss weight for xi]
                    bval = bval + g.gaussw(gi)*(Nxi*Nxi')*eljac(m,gi);
                    
                    %% Transformed Basis function gradients at xi
                    tg(:,3*(gi-1)+1:3*gi,m) = dNxi / Jac;
                end
                % Build up upper right part of symmetric mass matrix
                for j=1:eldofs
                    Mass(elem(j),elem(j:end)) = Mass(elem(j),elem(j:end)) + bval(j,j:end);
                end
            end
            this.Ngp = Ngpval;
            this.transgrad = tg;
            this.M = sparse(Mass + Mass' - diag(diag(Mass)));
            this.elem_detjac = eljac;
        end
                
        function plot(this, pm)
            if nargin < 2
                pm = PlotManager;%(false,1,2);
                pm.LeaveOpen = true;
            end
            
            h = pm.nextPlot('mass','Mass matrix','node','node');
            mesh(h,full(this.M));
            
            pm.nextPlot('mass_pattern','Mass matrix pattern','dof','dof');
            spy(this.M);
            
            if nargin < 2
                pm.done;
            end
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
                [pts, cubes] = geometry.Cube8Node.DemoGrid(ranges{k},ranges{k});
                g = geometry.Cube8Node(pts, cubes);
                tl = trilinear(g);
                tq = triquadratic(g.toCube20);
                res = res && norm(tl.elem_detjac-tq.elem_detjac,'inf') < 1e-14;
            end
        end
    end
    
end

