classdef Cube20Node < geometry.BaseGeometry
% Hexahedral geometry with 20 position nodes on each basic hexahedron/cube.
%
% Cube indexing on [-1, 1]³:
%
%           18---19---20 
%          / |       / |
%        16  |    17  |        Z+
%       /    11   /    12       |     Y+
%      13--14+---15    |        |    /
%      |     |   |     |        |   /
%      |     6---+7----8        |  /
%      |    /    |    /         | / 
%      9   4    10   5          |/
%      | /       | /            +---------X+
%      |/        |/
%      1----2----3 
%
% Corner indices: 1 3 6 8 13 15 18 20
%
%% Cube node positions:
% C1     -1    -1    -1
% E2      0    -1    -1
% C3      1    -1    -1
% E4     -1     0    -1
%       % 0     0    -1
% E5      1     0    -1
% C6     -1     1    -1
% E7      0     1    -1
% C8      1     1    -1
% E9     -1    -1     0 
%       % 0    -1     0
% E10     1    -1     0
%       %-1     0     0
%       % 0     0     0
%       % 1     0     0
% E11    -1     1     0
%       % 0     1     0
% E12     1     1     0
% C13    -1    -1     1
% E14     0    -1     1
% C15     1    -1     1
% E16    -1     0     1
%       % 0     0     1
% E17     1     0     1
% C18    -1     1     1
% E19     0     1     1
% C20     1     1     1
%
% without combinations 5,11,13,14,15,17 and 23 (they are neither on a
% corner or an edge due to 2 zero entries)

    methods
        
        function this = Cube20Node(pts, cubes)
            if nargin < 2
                [pts, cubes] = geometry.Cube20Node.DemoGrid;
            elseif size(unique(pts','rows'),1) ~= size(pts,2);
                error('Please provide unique points!');
            end
            this.Nodes = pts;
            this.Elements = cubes;
            this.DofsPerElement = 20;
            %EdgeIndices = [1 3 6 8 13 15 18 20];
            
            %% Compute edges
            e = int16.empty(0,2);
            for i=1:size(cubes,1)
                hlp = cubes(i,[1 2 1 4 1 9 2 3 3 5 3 10 4 6 5 8 6 7 ...
                    6 11 7 8 8 12 9 13 10 15 11 18 12 20 13 14 13 16 ...
                    14 15 15 17 16 18 17 20 18 19 19 20]);
                e(end+1:end+24,:) = reshape(hlp',2,[])';
            end
            e = unique(e,'rows','stable');
            this.Edges = e;
            
            %% Set Face indices
            this.MasterFaces = [  1 4 6 9 11 13 16 18
                            3 5 8 10 12 15 17 20
                            1:3 9 10 13:15
                            6:8 11 12 18:20
                            1:8
                            13:20];
             this.PatchFacesIdx = [ 1 4 6 11 18 16 13 9
                                    3 5 8 12 20 17 15 10
                                    1 2 3 10 15 14 13 9
                                    6 7 8 12 20 19 18 11
                                    1 2 3 5 8 7 6 4
                                    13 14 15 17 20 19 18 16];
            this.PatchesPerFace = 1;
            this.Faces = this.computeFaces;
            
            % Sanity checks
            this.ReverseAxesIndices = [3 2 1 5 4 8 7 6 10 9 12 11 15 14 13 17 16 20 19 18
                                       6:8 4 5 1:3 11 12 9 10 18:20 16 17 13:15
                                       13:20 9:12 1:8];
            this.OrientationCheckIndices(:,:,1) = [1:3; 6:8; 13:15; 18:20];
            this.OrientationCheckIndices(:,:,2) = [1 4 6; 3 5 8; 13 16 18; 15 17 20];
            this.OrientationCheckIndices(:,:,3) = [1 9 13; 3 10 15; 6 11 18; 8 12 20];
            this.checkOrientation;
        end
        
        function cube8 = toCube8Node(this)
            % Creates a 8 node cube geometry from this 20 node cube
            % geometry by simply leaving out the on-edge nodes.
            elems20 = this.Elements;
            nodes20 = this.Nodes;
            elems8 = elems20(:,[1 3 6 8 13 15 18 20]);
            usednodes = unique(elems8(:),'stable');
            nodes8 = nodes20(:,usednodes);
            invidx(usednodes) = 1:length(usednodes);
            elems8 = invidx(elems8);
            cube8 = geometry.Cube8Node(nodes8,elems8);
        end
        
    end
    
    methods(Static)
        function [pts, cubes] = DemoGrid(varargin)
            devperc = 0;
            if length(varargin) > 3
                devperc = varargin{4};
            end
            [pts, cubes] = geometry.Cube8Node.DemoGrid(varargin{1:min(length(varargin),3)});
            g8 = geometry.Cube8Node(pts, cubes);
            g20 = g8.toCube20Node;
            pts = g20.Nodes;
            cubes = g20.Elements;
            
            % Slightly deviate the grid
            if devperc > 0
                s = RandStream('mt19937ar','Seed',1);
                pts = pts + s.rand(size(pts))*devperc;
            end
        end
    end
    
    
end