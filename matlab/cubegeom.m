classdef cubegeom < handle
%

    properties
        % n x 3 position vector of nodes
        pts;
        
        % m x 8 index vector for all 8 points of m cubes
        cubes;
        
        % 2 x k index vector for edges between two points
        edges;
    end
    
    properties(Dependent)
        NumCubes;
        NumNodes;
        NumGaussp;
    end
    
    properties(SetAccess=private)
        gaussp;
        
        gaussw;
    end
    
    methods
        
        function this = cubegeom(pts, cubes)
            if nargin < 2
                [pts, cubes] = cubegeom.DemoCubeGrid;
            elseif size(unique(pts','rows'),1) ~= size(pts,2);
                error('Please provide unique points!');
            end
            this.pts = pts;
            this.cubes = cubes;

            % Init 27 Gauss points for 3-rule
            g = [-sqrt(3/5) 0 sqrt(3/5)];
            w = [5/9 8/9 5/9];
            [WX,WY,WZ] = meshgrid(w);
            [GX,GY,GZ] = meshgrid(g);
            W = WX.*WY.*WZ;
            this.gaussp = [GX(:) GY(:) GZ(:)]';
            this.gaussw = W(:);
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
            daspect([1 1 1]);
            if nargin < 2
                pm.done;
            end
        end
        
        function bounds = getBoundingBox(this, marginfac)
            if nargin < 2
                marginfac = 1;
            end
            m = min(this.pts,[],2)*marginfac;
            M = max(this.pts,[],2)*marginfac;
            bounds = [m(1) M(1) m(2) M(2) m(3) M(3)];
        end
    end
    
    methods
        function nc = get.NumCubes(this)
            nc = size(this.cubes,1);
        end
        
        function nc = get.NumNodes(this)
            nc = size(this.pts,2);
        end
        
        function nc = get.NumGaussp(this)
            nc = size(this.gaussp,2);
        end
        
    end
    
    methods(Static)
        function [pts, cubes] = DemoCubeGrid(xr,yr,zr)
            if nargin < 3
                zr = -1:1;
                if nargin < 2
                    yr = -1:1;
                    if nargin < 1
                        xr = -1:1;
                    end
                end
            end
            % Generate regular grid
            [X,Y,Z] = ndgrid(xr,yr,zr);
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
%             pts = pts + s.rand(size(pts))*.1;
        end
    end
    
    
end