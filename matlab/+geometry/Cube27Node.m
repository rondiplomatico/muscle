classdef Cube27Node < geometry.BaseGeometry
% Hexahedral geometry with 27 position nodes on each basic hexahedron/cube.
%
% Cube indexing on [-1, 1]³:
%
% Corner indices: 1 3 7 9 19 21 25 27
%
%% Cube node positions:
% C1     -1    -1    -1
% E2      0    -1    -1
% C3      1    -1    -1
% E4     -1     0    -1
% M5      0     0    -1
% E6      1     0    -1
% C7     -1     1    -1
% E8      0     1    -1
% C9      1     1    -1
% E10    -1    -1     0 
% M11     0    -1     0
% E12     1    -1     0
% M13    -1     0     0
% M14     0     0     0
% M15     1     0     0
% E16    -1     1     0
% M17     0     1     0
% E18     1     1     0
% C19    -1    -1     1
% E20     0    -1     1
% C21     1    -1     1
% E22    -1     0     1
% M23     0     0     1
% E24     1     0     1
% C25    -1     1     1
% E26     0     1     1
% C27     1     1     1

    methods
        
        function this = Cube27Node(nodes, elems)
            if nargin < 2
                [nodes, elems] = geometry.Cube27Node.DemoGrid;
            elseif size(unique(nodes','rows'),1) ~= size(nodes,2);
                error('Please provide unique points!');
            end
            this.Nodes = nodes;
            this.Elements = elems;
            this.DofsPerElement = 27;
            %EdgeIndices = [1 3 7 9 19 21 25 27];
            
            %% Compute edges
            e = int16.empty(0,2);
            for i=1:size(elems,1)
                hlp = elems(i,[1 2 4 5 7 8 2 3 5 6 8 9 1 4 2 5 3 6 4 7 5 8 6 9 ... % bottom
                    1 10 2 11 3 12 10 19 11 20 12 21 10 11 11 12 19 20 20 21 ... % front 
                    12 15 15 18 21 24 24 27 6 15 15 24 9 18 18 27 ... % right side
                    4 13 13 22 7 16 16 25 10 13 13 16 19 22 22 25 ... % left side
                    8 17 17 26 16 17 17 18 25 26 26 27 ... % back side
                    20 23 23 26 22 23 23 24 ... % top
                    ]); % inner 5 14 13 14 11 14 14 15 14 17 14 23
                e(end+1:end+48,:) = reshape(hlp',2,[])';
            end
            e = unique(e,'rows','stable');
            this.Edges = e;
            
            %% Set Face indices
            this.MasterFaces = [1:3:25
                            3:3:27
                            1:3 10:12 19:21
                            7:9 16:18 25:27
                            1:9
                            19:27];
            % Face node indices for patch objects
            this.PatchFacesIdx = [1 4 13 10; 4 7 16 13; 13 16 25 22; 10 13 22 19;
                3 6 15 12; 6 9 18 15; 12 15 24 21; 15 18 27 24;
                1 2 11 10; 2 3 12 11; 10 11 20 19; 11 12 21 20;
                7 8 17 16; 8 9 18 17; 16 17 26 25; 17 18 27 26;
                1 2 5 4; 4 5 8 7; 2 3 6 5; 5 6 9 8;
                19 20 23 22; 20 21 24 23; 22 23 26 25; 23 24 27 26];
            this.PatchesPerFace = 4;
            this.Faces = this.computeFaces;
            
            this.ReverseAxesIndices = [3 2 1 6 5 4 9 8 7 12 11 10 15 14 13 18 17 16 21 20 19 24 23 22 27 26 25
                                       7:9 4:6 1:3 16:18 13:15 10:12 25:27 22:24 19:21
                                       19:27 10:18 1:9];
            
            this.OrientationCheckIndices(:,:,1) = [1:3; 4:6; 7:9; 10:12; 13:15; 16:18; 19:21; 22:24; 25:27];
            this.OrientationCheckIndices(:,:,2) = [1 4 7; 2 5 8; 3 6 9; 10 13 16; 11 14 17;12 15 18; 19 22 25; 20 23 26; 21 24 27];
            this.OrientationCheckIndices(:,:,3) = [1 10 19; 2 11 20; 3 12 21; 4 13 22; 5 14 23; 6 15 24; 7 16 25; 8 17 26; 9 18 27];
            this.checkOrientation;
        end
        
        function cube8 = toCube8Node(this)
            % Creates a 8 node cube geometry from this 20 node cube
            % geometry by simply leaving out the on-edge nodes.
            elems8 = this.Elements(:,[1 3 7 9 19 21 25 27]);
            usednodes = unique(elems8(:),'stable');
            invidx(usednodes) = 1:length(usednodes);
            elems8 = invidx(elems8);
            cube8 = geometry.Cube8Node(this.Nodes(:,usednodes),elems8);
        end
        
        function cube20 = toCube20Node(this)
            elems20 = this.Elements(:,[1:4 6:10 12 16 18:22 24:27]);
            usednodes = unique(elems20(:),'stable');
            invidx(usednodes) = 1:length(usednodes);
            elems20 = invidx(elems20);
            cube20 = geometry.Cube20Node(this.Nodes(:,usednodes),elems20);
        end
        
        function swapYZ(this)
            n = this.Nodes;
            n([2 3],:) = n([3 2],:);
            this.Nodes = n;
            this.Elements = this.Elements(:,[1:3 10:12 19:21 4:6 13:15 22:24 7:9 16:18 25:27]);
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
            g27 = g8.toCube27Node;
            pts = g27.Nodes;
            cubes = g27.Elements;
            
            % Slightly deviate the grid
            if devperc > 0
                s = RandStream('mt19937ar','Seed',1);
                pts = pts + s.rand(size(pts))*devperc;
            end
        end
    end
    
    
end