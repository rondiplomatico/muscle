classdef BaseGeometry < handle
%

    properties(SetAccess=protected)
        % n x 3 position vector of nodes
        Nodes;
        
        % m x p index vector for all p nodes of m elements
        Elements;
        
        % 2 x k index vector for edges between two points
        Edges;
        
        DofsPerElement;
        
        %CornerIndices;
    end
    
    properties(Dependent)
        NumElements;
        NumNodes;
        GaussPointsPerElem;
    end
    
    properties(SetAccess=private)
        gaussp;
        
        gaussw;
    end
    
    methods
        
        function this = BaseGeometry
            % Init 27 Gauss points for 3-rule
            g = [-sqrt(3/5) 0 sqrt(3/5)];
            w = [5/9 8/9 5/9];
            [WX,WY,WZ] = meshgrid(w);
            [GX,GY,GZ] = meshgrid(g);
            W = WX.*WY.*WZ;
            this.gaussp = [GX(:) GY(:) GZ(:)]';
            this.gaussw = W(:);
            
%             %% Compute node to cube indices
%             pv = cell(np,1);
%             for i=1:np
%                 pv{i} = find(sum(el == i,2));
%             end
%             this.pts_cubes = pv;
        end
        
        function plot(this, withdofnr, pm)
            if nargin < 3
                if nargin < 2
                    withdofnr = false;
                end
                pm = PlotManager;%(false,1,2);
                pm.LeaveOpen = true;
            end
            
            p = this.Nodes;
            h = pm.nextPlot('geometry',sprintf('Geometry (%s)',class(this)), 'x [mm]', 'y [mm]');
            plot3(h,p(1,:),p(2,:),p(3,:),'k.','MarkerSize',14);
            hold(h,'on');
            for k = 1:this.NumElements
                el = this.Elements(k,:);
                center = sum(p(:,el),2)/this.DofsPerElement;
                text(center(1),center(2),center(3),sprintf('#%d',k),'Parent',h,'Color','m');
                % Plot local numbering for first element
                if k == 1
                    off = .04;
                    for i = 1:this.DofsPerElement
                        text(off+p(1,el(i)),off+p(2,el(i)),off+p(3,el(i)),sprintf('%d',i),'Parent',h,'Color','r');
                    end
                end
            end
            eg = this.Edges;
            for i = 1:size(eg,1)
%                 if any(el == eg(i,1)) && any(el == eg(i,2))
                    plot3(h,[p(1,eg(i,1)) p(1,eg(i,2))],[p(2,eg(i,1)) p(2,eg(i,2))],[p(3,eg(i,1)) p(3,eg(i,2))],'r');
%                 end
            end
            if withdofnr
                for k = 1:this.NumNodes
                    text(p(1,k),p(2,k),p(3,k),sprintf('%d',k),'Parent',h,'Color','k');
                end
            end
            zlabel('z [mm]');
            view(h,[52 30]);
            daspect([1 1 1]);
                        
            if nargin < 2
                pm.done;
            end
        end
        
        function bounds = getBoundingBox(this, marginfac)
            if nargin < 2
                marginfac = 1;
            end
            m = min(this.Nodes,[],2)*marginfac;
            M = max(this.Nodes,[],2)*marginfac;
            bounds = [m(1) M(1) m(2) M(2) m(3) M(3)];
        end
        
        function commonidx = getCommonNodesWith(this, other)
            commonidx = Utils.findVecInMatrix(this.Nodes,other.Nodes);
        end
           
    end
    
    methods
        function nc = get.NumElements(this)
            nc = size(this.Elements,1);
        end
        
        function nc = get.NumNodes(this)
            nc = size(this.Nodes,2);
        end
        
        function nc = get.GaussPointsPerElem(this)
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