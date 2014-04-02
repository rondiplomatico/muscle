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
            
            %% Mass matrix
            % Iterate each cube
            gp = g.gaussp;
            %             gjac = zeros(size(c,1),size(gp,2));
            Mass = zeros(np,np);
            eldofs = size(el,2);
            for m = 1:size(el,1)
                bval = zeros(eldofs,eldofs);
                elem = el(m,:);
                for gi = 1:size(gp,2)
                    xi = gp(:,gi); % gauss point \xi
                    J = det(dof(:,elem)*this.gradN(xi));
                    %                     gjac(m,g) = det(J);
                    bval = bval + g.gaussw(gi)*this.N(xi)*this.N(xi)'*J;
                end
                for j=1:eldofs
                    Mass(elem(j),elem(j:end)) = Mass(elem(j),elem(j:end)) + bval(j,j:end);
                end
            end
            this.M = sparse(Mass + Mass');
            %             this.cube_detjac = gjac;
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
    
end

