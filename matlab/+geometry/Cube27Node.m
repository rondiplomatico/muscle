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
        
        function this = Cube27Node(pts, cubes)
            if nargin < 2
                [pts, cubes] = geometry.Cube27Node.DemoGrid;
            elseif size(unique(pts','rows'),1) ~= size(pts,2);
                error('Please provide unique points!');
            end
            this.Nodes = pts;
            this.Elements = cubes;
            this.DofsPerElement = 27;
            %EdgeIndices = [1 3 7 9 19 21 25 27];
            
            %% Compute edges
            e = int16.empty(0,2);
            for i=1:size(cubes,1)
                hlp = cubes(i,[1 2 4 5 7 8 2 3 5 6 8 9 1 4 2 5 3 6 4 7 5 8 6 9 ... % bottom
                    1 10 2 11 3 12 10 19 11 20 12 21 10 11 11 12 19 20 20 21 ... % front 
                    12 15 15 18 21 24 24 27 6 15 15 24 9 18 18 27 ... % right side
                    4 13 13 22 7 16 16 25 10 13 13 16 19 22 22 25 ... % left side
                    8 17 17 26 16 17 17 18 25 26 26 27 ... % back side
                    20 23 23 26 22 23 23 24 ... % top
                    5 14 13 14 11 14 14 15 14 17 14 23]); % inner
                e(end+1:end+54,:) = reshape(hlp',2,[])';
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
            this.Faces = this.computeFaces;
        end
        
        function cube8 = toCube8Node(this)
            % Creates a 8 node cube geometry from this 20 node cube
            % geometry by simply leaving out the on-edge nodes.
            elems20 = this.Elements;
            nodes20 = this.Nodes;
            elems8 = elems20(:,[1 3 7 9 19 21 25 27]);
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