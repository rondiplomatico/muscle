classdef cubegeom < handle
%

%% Cube indexing:
%  /7---8 1: (-1,-1,-1)
% 3-+-4/| 2: ( 1,-1,-1)
% | 5-+-6 3: (-1, 1,-1) 
% 1---2/  4: ( 1, 1,-1)
%         5: (-1,-1, 1)
%         6: ( 1,-1, 1)
%         7: (-1, 1, 1)
%         8: ( 1, 1, 1)
    
    properties
        % n x 3 position vector of nodes
        pts;
        
        % m x 8 index vector for all 8 points of m cubes
        cubes;
        
        % 2 x k index vector for edges between two points
        edges;
        
        % cell array of dim n containing indices of all cubes that are
        % adjacent to the n-th point
        pts_cubes;
        
        % m x g The determinants of the deformation jacobians at the gauss
        % points
%         cube_detjac;
        
        % n x 1 cell array containing the indices of the neighboring points
        % of the n-th point
%         neighbors;
        
        % The mass matrix
        M;
    end
    
    properties(Dependent)
        NumCubes;
        NumNodes;
    end
    
    properties(SetAccess=private)
        N = @(x)[   (1-x(1,:)).*(1-x(2,:)).*(1-x(3,:));... % 1
            (1+x(1,:)).*(1-x(2,:)).*(1-x(3,:));...
            (1-x(1,:)).*(1+x(2,:)).*(1-x(3,:));... % 3
            (1+x(1,:)).*(1+x(2,:)).*(1-x(3,:));...
            (1-x(1,:)).*(1-x(2,:)).*(1+x(3,:));... % 5
            (1+x(1,:)).*(1-x(2,:)).*(1+x(3,:));...
            (1-x(1,:)).*(1+x(2,:)).*(1+x(3,:));... % 7
            (1+x(1,:)).*(1+x(2,:)).*(1+x(3,:));]/8;

        gradN = @(x)[-(1-x(2,:)).*(1-x(3,:)) -(1-x(1,:)).*(1-x(3,:)) -(1-x(1,:)).*(1-x(2,:));... % 1
            (1-x(2,:)).*(1-x(3,:)) -(1+x(1,:)).*(1-x(3,:)) -(1+x(1,:)).*(1-x(2,:));...
            -(1+x(2,:)).*(1-x(3,:)) (1-x(1,:)).*(1-x(3,:)) -(1-x(1,:)).*(1+x(2,:));... % 3
            (1+x(2,:)).*(1-x(3,:))  (1+x(1,:)).*(1-x(3,:)) -(1+x(1,:)).*(1+x(2,:));...
            -(1-x(2,:)).*(1+x(3,:)) -(1-x(1,:)).*(1+x(3,:)) (1-x(1,:)).*(1-x(2,:));... % 5
            (1-x(2,:)).*(1+x(3,:)) -(1+x(1,:)).*(1+x(3,:)) (1+x(1,:)).*(1-x(2,:));...
            -(1+x(2,:)).*(1+x(3,:)) (1-x(1,:)).*(1+x(3,:)) (1-x(1,:)).*(1+x(2,:));... % 7
            (1+x(2,:)).*(1+x(3,:))  (1+x(1,:)).*(1+x(3,:)) (1+x(1,:)).*(1+x(2,:))]/8;
        
        gaussp;
        
        gaussw;
    end
    
    methods
        
        function this = cubegeom(pts, cubes)
            if nargin < 2
                [pts, cubes] = cubegeom.DemoCubeGrid;
            end
            this.pts = pts;
            this.cubes = cubes;
            s = RandStream('mt19937ar','Seed',1);
%             this.gaussp = 2*(s.rand(3,5)-.5);
%             this.gaussw = s.rand(1,5);
            this.gaussp = [0 0 0]';
            this.gaussw = 1;
            
            this.precomputeStuff;
        end
        
        function precomputeStuff(this)
            %% Compute edges from cubes.
            % Mind the corner numbering of the cubes! (see above)
            e = int16.empty(0,2);
            c = this.cubes;
            p = this.pts;
            np = size(p,2);
            for i=1:size(c,1)
                hlp = c(i,[1 2 1 3 1 5 3 4 2 4 4 8 3 7 ...
                           8 7 5 7 6 2 6 5 6 8]);
                e(end+1:end+12,:) = reshape(hlp',2,[])';
            end
            e = unique(e,'rows');
            this.edges = e;
            
            %% Compute neighbors
%             n = cell(np,1);
%             for i=1:np
%                 n{i} = unique([e(e(:,1) == i,2); e(e(:,2) == i,1)])';
%             end
%             this.neighbors = n;
            
            %% Compute node to cube indices
            pv = cell(np,1);
            for i=1:np
               pv{i} = find(sum(c == i,2));
            end
            this.pts_cubes = pv;
            
            %% Mass matrix
            % Iterate each cube
            gp = this.gaussp;
%             gjac = zeros(size(c,1),size(gp,2));
            Mass = zeros(np,np);
            for m = 1:size(c,1)
                bval = zeros(8,8);
                cube = c(m,:);
                for g = 1:size(gp,2)
                    xi = gp(:,g); % gauss point \xi
                    J = det(p(:,cube)*this.gradN(xi));
%                     gjac(m,g) = det(J);
                    bval = bval + this.gaussw(g)*this.N(xi)*this.N(xi)'*J;
                end
                for j=1:8
                    Mass(cube(j),cube(j:end)) = Mass(cube(j),cube(j:end)) + bval(j,j:end);
                end
            end
            this.M = sparse(Mass + Mass');
%             this.cube_detjac = gjac;
        end
        
        function plot(this, pm)
            if nargin < 2
                pm = PlotManager(false,1,2);
                pm.LeaveOpen = true;
            end
            p = this.pts;
            h = pm.nextPlot('grid','CubeGrid','x','y');
            plot3(h,p(1,:),p(2,:),p(3,:),'k.','MarkerSize',14);
            e = this.edges;
            hold(h,'on');
            for k=1:size(e,1)
                plot3(h,p(1,[e(k,1) e(k,2)]),p(2,[e(k,1) e(k,2)]),p(3,[e(k,1) e(k,2)]),'r');
            end
            h = pm.nextPlot('mass','Mass matrix','node','node');
            mesh(h,full(this.M));
            
            if nargin < 2
                pm.done;
            end
        end
    end
    
    methods
        function nc = get.NumCubes(this)
            nc = size(this.cubes,1);
        end
        
        function nc = get.NumNodes(this)
            nc = size(this.pts,2);
        end
    end
    
    methods(Static)
%         function [pts, edges] = DemoCubeGrid
%             % Generate regular grid
%             [X,Y,Z] = ndgrid(-1:1,-1:1,-1:1);
%             pts = [X(:) Y(:) Z(:)];
%             
%             % Hack: find all points with distance one
%             sq = sum(pts.*pts,2)';
%             n = size(pts,1);
%             r = ((ones(n,1)*sq)' + ones(n,1)*sq - 2*(pts*pts'));
%             [i,j] = find(tril(r) == 1);
%             edges = [j i];
%             
%             % Slightly deviate the grid
%             pts = pts + rand(size(pts))*.05;
%         end

        function res = test_CubeGeom
            c = cubegeom;
            [X,Y,Z] = ndgrid(-1:2:1,-1:2:1,-1:2:1);
            p = [X(:) Y(:) Z(:)]';
            if ~isequal(c.N(p),eye(8))
                error('Basis function mismatch');
            end
            res = true;
        end
        
        function [pts, cubes] = DemoCubeGrid
            % Generate regular grid
            [X,Y,Z] = ndgrid(-1:1,-2:1,-1:3);
            pts = [X(:) Y(:) Z(:)]';
            
            cubes = double.empty(0,8);
            for i = 1:size(X,1)-1
                for j = 1:size(X,2)-1
                    for k = 1:size(X,3)-1
                        hx = X([i i+1],[j j+1],[k k+1]);
                        hy = Y([i i+1],[j j+1],[k k+1]);
                        hz = Z([i i+1],[j j+1],[k k+1]);
                        cubes(end+1,:) = Utils.findVecInMatrix(pts,[hx(:) hy(:) hz(:)]');%#ok
                    end
                end
            end
            
            % Slightly deviate the grid
%             s = RandStream('mt19937ar','Seed',1);
%             pts = pts + s.rand(size(pts))*.2;
        end
    end
    
    
end