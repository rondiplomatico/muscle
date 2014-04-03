classdef trilinear < fembase
    %TRILINEAR Summary of this class goes here
    %   Detailed explanation goes here
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
        
        % cell array of dim n containing indices of all cubes that are
        % adjacent to the n-th point
        pts_cubes;
    end
    
    methods
        function this = trilinear(geo)
            if nargin < 1
                geo = cubegeom;
            end
            
            this = this@fembase(geo);
            
            this.init;
        end
        
        function init(this)
            % Trilinear basis functions
            this.N = @(x)[   (1-x(1,:)).*(1-x(2,:)).*(1-x(3,:));... % 1
                (1+x(1,:)).*(1-x(2,:)).*(1-x(3,:));...
                (1-x(1,:)).*(1+x(2,:)).*(1-x(3,:));... % 3
                (1+x(1,:)).*(1+x(2,:)).*(1-x(3,:));...
                (1-x(1,:)).*(1-x(2,:)).*(1+x(3,:));... % 5
                (1+x(1,:)).*(1-x(2,:)).*(1+x(3,:));...
                (1-x(1,:)).*(1+x(2,:)).*(1+x(3,:));... % 7
                (1+x(1,:)).*(1+x(2,:)).*(1+x(3,:));]/8;
            
            this.gradN = @(x)[-(1-x(2,:)).*(1-x(3,:)) -(1-x(1,:)).*(1-x(3,:)) -(1-x(1,:)).*(1-x(2,:));... % 1
                (1-x(2,:)).*(1-x(3,:)) -(1+x(1,:)).*(1-x(3,:)) -(1+x(1,:)).*(1-x(2,:));...
                -(1+x(2,:)).*(1-x(3,:)) (1-x(1,:)).*(1-x(3,:)) -(1-x(1,:)).*(1+x(2,:));... % 3
                (1+x(2,:)).*(1-x(3,:))  (1+x(1,:)).*(1-x(3,:)) -(1+x(1,:)).*(1+x(2,:));...
                -(1-x(2,:)).*(1+x(3,:)) -(1-x(1,:)).*(1+x(3,:)) (1-x(1,:)).*(1-x(2,:));... % 5
                (1-x(2,:)).*(1+x(3,:)) -(1+x(1,:)).*(1+x(3,:)) (1+x(1,:)).*(1-x(2,:));...
                -(1+x(2,:)).*(1+x(3,:)) (1-x(1,:)).*(1+x(3,:)) (1-x(1,:)).*(1+x(2,:));... % 7
                (1+x(2,:)).*(1+x(3,:))  (1+x(1,:)).*(1+x(3,:)) (1+x(1,:)).*(1+x(2,:))]/8;
            
            this.nodes = this.geo.pts;
            this.elems = this.geo.cubes;
            
            init@fembase(this);
            
            %% Compute edges
            e = int16.empty(0,2);
            for i=1:size(this.elems,1)
                hlp = this.elems(i,[1 2 1 3 1 5 3 4 2 4 4 8 3 7 ...
                    8 7 5 7 6 2 6 5 6 8]);
                e(end+1:end+12,:) = reshape(hlp',2,[])';
            end
            e = unique(e,'rows');
            this.edges = e;
        end
        
        
    end
    
    methods(Static)
        function res = test_TrilinearBasisFun
            q = trilinear;
            res = fembase.test_BasisFun(q);
            
            % test for correct basis function values on nodes
            [X,Y,Z] = ndgrid(-1:2:1,-1:2:1,-1:2:1);
            p = [X(:) Y(:) Z(:)]';
            res = res && isequal(q.N(p),eye(8));
        end
    end
    
end

