classdef Cube20Node < geometry.BaseGeometry
%% Cube indexing:
%
%           18---19---20 
%          / |       / |
%        /   |     /   |        Y+
%       /    16   /    17       |     Z+
%      6----7+---8     |        |    /
%      |     |   |     |        |   /
%      |     13--+14---15       |  /
%      |    /    |    /         | / 
%      4   9     5  10          |/
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
                [pts, cubes] = Cube20Node.DemoGrid;
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
            e = unique(e,'rows');
            this.Edges = e;
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
            [pts, cubes] = geometry.Cube8Node.DemoGrid(varargin{:});
            
            % Slightly deviate the grid
%             s = RandStream('mt19937ar','Seed',1);
%             pts = pts + s.rand(size(pts))*.1;

            g8 = geometry.Cube8Node(pts, cubes);
            g20 = g8.toCube20Node;
            pts = g20.Nodes;
            cubes = g20.Elements;
        end
    end
    
    
end