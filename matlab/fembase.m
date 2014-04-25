classdef fembase < handle
    %FEMBASE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        geo;
        
        N;
        
        % Values of basis functions on gauss points
        Ngp;
        
        gradN;
        
        % The mass matrix
        M;
        
        nodes;
        
        elems;
        
        elem_detjac;
        
        edges;
        
        % transformed basis function gradients, stored in eldofs x gp*3
        % x element matrix (accessible in 20x3-chunks for each gauss point and element)
        transgrad;
    end
    
    properties(Dependent)
        NumNodes;
        NumElems;
        DofsPerElement;
    end
    
    properties(SetAccess=protected)
        % The indices of edges
        EdgeIndices = [1 3 6 8 13 15 18 20];
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
            dof = this.nodes;
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
                    Jac = this.nodes(:,elem)*dNxi;
                    
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
            this.M = sparse(Mass + Mass');
            this.elem_detjac = eljac;
        end
        
%         function J = getJacobian(this, elem, xi)
%             J = this.nodes(:,elem)*this.gradN(xi);
%         end
        
        function plot(this, pm)
            if nargin < 2
                pm = PlotManager;%(false,1,2);
                pm.LeaveOpen = true;
            end
            this.geo.plot(pm);
            
            p = this.nodes;
            h = pm.nextPlot('nodes','Grid nodes','x','y');
            plot3(h,p(1,:),p(2,:),p(3,:),'k.','MarkerSize',14);
            for k = 1:this.NumElems
                el = this.elems(k,:);
                center = sum(p(:,el),2)/this.DofsPerElement;
                text(center(1),center(2),center(3),sprintf('#%d',k),'Parent',h,'Color','m');
                % Plot local numbering for first element
                if k == 1
                    off = .04;
                    for i = 1:this.DofsPerElement
                        text(off+p(1,el(i)),off+p(2,el(i)),off+p(3,el(i)),sprintf('%d',i),'Parent',h,'Color','r');
                    end
                    eg = this.edges;
                    hold(h,'on');
                    for i = 1:size(eg,1)
                        if any(el == eg(i,1)) && any(el == eg(i,2))
                            plot3(h,[p(1,eg(i,1)) p(1,eg(i,2))],[p(2,eg(i,1)) p(2,eg(i,1))],[p(3,eg(i,1)) p(3,eg(i,2))],'r');
                        end
                    end
                end
            end
            for k = 1:this.NumNodes
                text(p(1,k),p(2,k),p(3,k),sprintf('%d',k),'Parent',h,'Color','k');
            end
            view(h,[52 30]);
            
%             h = pm.nextPlot('mass','Mass matrix','node','node');
%             mesh(h,full(this.M));
%             
%             pm.nextPlot('mass_pattern','Mass matrix pattern','dof','dof');
%             spy(this.M);
            
            if nargin < 2
                pm.done;
            end
        end
    end
    
    methods
        function value = get.NumNodes(this)
            value = size(this.nodes,2);
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
                res = res && norm(tl.elem_detjac-tq.elem_detjac,'inf') < 1e-14;
            end
        end
    end
    
end

