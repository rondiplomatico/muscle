classdef Cube8Node < geometry.BaseGeometry
%% Cube indexing:
%  /7---8 1: (-1,-1,-1)
% 5-+-6/| 2: ( 1,-1,-1)
% | 3-+-4 3: (-1, 1,-1) 
% 1---2/  4: ( 1, 1,-1)
%         5: (-1,-1, 1)
%         6: ( 1,-1, 1)
%         7: (-1, 1, 1)
%         8: ( 1, 1, 1)

    methods
        
        function this = Cube8Node(pts, cubes)
            if nargin < 2
                [pts, cubes] = geometry.Cube8Node.DemoGrid;
            elseif size(unique(pts','rows'),1) ~= size(pts,2);
                error('Please provide unique points!');
            end
            this.Nodes = pts;
            this.Elements = cubes;
            this.DofsPerElement = 8;
            
            %% Compute edges
            e = int16.empty(0,2);
            for i=1:size(cubes,1)
                hlp = cubes(i,[1 2 1 3 1 5 3 4 2 4 4 8 3 7 ...
                    8 7 5 7 6 2 6 5 6 8]);
                e(end+1:end+12,:) = reshape(hlp',2,[])';
            end
            e = unique(e,'rows','stable');
            this.Edges = e;
            
            %% Set Face indices
            this.MasterFaces = [1 3 5 7
                            2 4 6 8
                            1 2 5 6
                            3 4 7 8
                            1 2 3 4
                            5 6 7 8];
            this.PatchFacesIdx = [1 3 7 5
                            2 4 8 6
                            1 2 6 5
                            3 4 8 7
                            1 2 4 3
                            5 6 8 7];
            this.PatchesPerFace = 1;
            this.Faces = this.computeFaces;
        end
        
        function cube20 = toCube20Node(this)
            % Creates a 20 node cube geometry from this 8 node cube
            % geometry by linear interpolation between corners.
            elems8 = this.Elements;
            nodes8 = this.Nodes;
            nc = size(elems8,1);
            % Transformation matrix for corners to corner+edges locations
            i = [1 2  2  3  4  4  5  5 6  7  7 8  9  9 10 10 11 11 12 12 13 14 14 15 16 16 17 17 18 19 19 20];  
            j = [1 1  2  2  1  3  2  4 3  3  4 4  1  5 2  6  3  7  4  8  5  5  6  6  5  7  6  8  7  7  8  8];
            s = [1 .5 .5 1 .5 .5 .5 .5 1 .5 .5 1 .5 .5 .5 .5 .5 .5 .5 .5 1  .5 .5 1  .5 .5 .5 .5 1  .5 .5 1];
            T = sparse(i,j,s,20,8);
            nodes20 = zeros(3,nc*20);

            % Iterate all cubes and collect nodes
            for cidx = 1:nc
                pos = 20*(cidx-1)+1:20*cidx;
                cube = elems8(cidx,:);
                nodes20(:,pos) = nodes8(:,cube)*T';
            end
            [nodes20, ~, elems20] = unique(nodes20','rows','stable');
            cube20 = geometry.Cube20Node(nodes20',reshape(elems20,20,[])');
        end
        
        function cube27 = toCube27Node(this)
            % Creates a 20 node cube geometry from this 8 node cube
            % geometry by linear interpolation between corners.
            elems8 = this.Elements;
            nodes8 = this.Nodes;
            nc = size(elems8,1);
            h = .5;
            v = .25;
            % Transformation matrix for corners to corner+edges locations
            % Nodes 1-13
            i = [1 2  2  3  4  4  5 5 5 5 6 6 7 8 8 9 10 10 11 11 11 11 12 12 13 13 13 13];
            j = [1 1  2  2  1  3  1 2 3 4 2 4 3 3 4 4 1  5  1  2  5  6  2  6  1  3  5  7]; 
            s = [1 h  h  1  h  h  v v v v h h 1 h h 1 .5 .5 v  v  v  v  .5 .5 v  v  v  v];
            % Node 14
            i = [i ones(1,8)*14];
            j = [j 1:8];
            s = [s ones(1,8)/8];
            % Nodes 15-27
            i = [i 15 15 15 15 16 16 17 17 17 17 18 18 19 20 20 21 22 22 23 23 23 23 24 24 25 26 26 27];
            j = [j 2  4  6  8  3  7  3  4  7  8  4  8  5  5  6  6  5  7  5  6  7  8  6  8  7  7  8  8]; 
            s = [s v  v  v  v  h  h  v  v  v  v  h  h  1  h  h  1  .5 .5 v  v  v  v  .5 .5 1  h  h  1];
            T = sparse(i,j,s,27,8);
            nodes27 = zeros(3,nc*27);

            % Iterate all cubes and collect nodes
            for cidx = 1:nc
                pos = 27*(cidx-1)+1:27*cidx;
                cube = elems8(cidx,:);
                nodes27(:,pos) = nodes8(:,cube)*T';
            end
            [nodes27, ~, elems27] = unique(nodes27','rows','stable');
            cube27 = geometry.Cube27Node(nodes27',reshape(elems27,27,[])');
        end

    end
    
    methods(Static)
        function [pts, cubes] = DemoGrid(xr,yr,zr,devperc)
            if nargin < 4
                devperc = 0;
                if nargin < 3
                    zr = [-1 1];
                    if nargin < 2
                        yr = [-1 1];
                        if nargin < 1
                            xr = [-1 1];
                        end
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
            if devperc > 0
                s = RandStream('mt19937ar','Seed',1);
                pts = pts + s.rand(size(pts))*devperc;
            end
        end
    end
    
    
end