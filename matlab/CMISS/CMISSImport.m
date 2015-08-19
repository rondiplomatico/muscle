classdef CMISSImport < handle
    %IMPORT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        Dir = fileparts(which(mfilename));
    end
    
    properties
        % Flag that indicates if the nodes should be centered around 0
        % after loading
        CenterNodes = true;
    end
    
    methods(Static)
        function [geo8, geo20, geo27] = readEntireTA
            ci = CMISSImport;
            name = 'EntireTA';
            [nodes, nodeidx] = ci.readNodes(name);
            
            elems27(1,:) = [6,21,103,31,33,155,8,22,104,46,47,182,48,49,184,70,71,230,13,23,125,32,34,156,15,24,127];
            elems27(2,:) = [13,23,125,32,34,156,15,24,127,54,55,202,56,57,204,74,75,240,9,25,121,35,37,169,11,26,123];
            elems27(3,:) = [9,25,121,35,37,169,11,26,123,50,51,198,52,53,200,72,73,238,17,27,129,36,38,170,19,28,131];
            elems27(4,:) = [17,27,129,36,38,170,19,28,131,58,59,206,60,61,208,76,77,242,120,194,101,166,196,165,117,236,102];
            elems27(5,:) = [103,137,105,155,157,173,104,138,106,182,183,214,184,185,215,230,231,246,125,139,126,156,158,174,127,140,128];
            elems27(6,:) = [109,141,113,159,161,175,110,142,114,186,187,216,188,189,217,232,233,247,133,143,134,160,162,176,135,144,136];
            elems27(7,:) = [111,145,115,163,164,177,112,146,116,190,191,218,192,193,219,234,235,248,103,137,105,155,157,173,104,138,106];
            elems27(8,:) = [101,147,107,165,167,178,102,148,108,194,195,220,196,197,221,236,237,249,120,149,119,166,168,179,117,150,118];
            elems27(9,:) = [121,151,122,169,171,180,123,152,124,198,199,222,200,201,223,238,239,250,129,153,130,170,172,181,131,154,132];
            elems27(10,:) = [125,139,126,156,158,174,127,140,128,202,203,224,204,205,225,240,241,251,121,151,122,169,171,180,123,152,124];
            elems27(11,:) = [129,153,130,170,172,181,131,154,132,206,207,226,208,209,227,242,243,252,101,147,107,165,167,178,102,148,108];
            elems27(12,:) = [133,143,134,160,162,176,135,144,136,210,211,228,212,213,229,244,245,253,111,145,115,163,164,177,112,146,116];
            
            [nodes27, elems27] = ci.merge(nodes, nodeidx, elems27);
            geo27 = fem.geometry.Cube27Node(nodes27, elems27);
            geo27.swapYZ;
            
            geo8 = geo27.toCube8Node;
            geo20 = geo27.toCube20Node;
            
            save(fullfile(CMISSImport.Dir,name),'geo8','geo20','geo27');
        end
    end
    
    methods
        
        function [geo8, geo20, geo27] = import(this, name)
            [nodes, nodeidx] = this.readNodes(name);
            [elems8, elems20, elems27] = this.readElems(name);
            [nodes8, elems8] = this.merge(nodes, nodeidx, elems8);
            [nodes20, elems20] = this.merge(nodes, nodeidx, elems20);
            [nodes27, elems27] = this.merge(nodes, nodeidx, elems27);
            geo8 = fem.geometry.Cube8Node(nodes8, elems8);
            geo20 = fem.geometry.Cube20Node(nodes20, elems20);
            geo27 = fem.geometry.Cube27Node(nodes27, elems27);
            save(fullfile(CMISSImport.Dir,name),'geo8','geo20','geo27');
        end
        
        function [nodes, nodeidx] = readNodes(this, name)
            filename = fullfile(CMISSImport.Dir,[name '.ipnode']);
            [nodes, nodenr] = this.rawReadNodes(filename);
            empt = cellfun(@(x)isempty(x),nodenr);
            nodes(empt) = [];
            nodenr(empt) = [];
            
            nodenr(1:9) = [];
            nodeidx = eval(['[' sprintf('%s ',nodenr{1:4:end}) ']']);
            nodes(1:10) = [];
            nodes(4:4:end) = [];
            nodes = reshape(nodes,3,[]);
            if this.CenterNodes
                nodes = nodes - repmat(mean(nodes,2),1,size(nodes,2));
            end
        end
        
        function [elems8, elems20, elems27] = readElems(this, name)
            filename = fullfile(CMISSImport.Dir,[name '.ipelem']);
            data = this.rawReadElems(filename);
            data(1:5,:) = [];
            %elemidx = data(1:9:end,5);
            
            %% Read 27 points for quadratic case
            
            elems27 = [data(6:9:end,9:end) data(7:9:end,1:11)];
            % Modified read for ip data from M. Sprenger
%             elems27 = [data(6:9:end,9:end-1) data(7:9:end,1:12)];
            
            % Remove the inner/face points
            elems20 = elems27;
            elems20(:,[5 11 13:15 17 23]) = [];
            
            %% Read 8 points for linear case
            elems8 = data(8:9:end,9:16);
        end
        
        function [nodes, elems] = merge(~, nodes, nodeidx, elems)
            %% Determine effectively used nodes
            effidx = unique(elems(:));
            [~, idx] = intersect(nodeidx, effidx);
            % Remove unused
            nodes = nodes(:,idx);
            % Fit element node reference indices
            invidx(effidx) = 1:length(effidx);
            elems = invidx(elems);
            
            % Sort the node indices per element in the order of z
            % increasing "outer", y increasing "mid" and x increasing "inner"
%             for m=1:size(elems,1)
%                 idx = sortXYZ(nodes(:,elems(m,:)));
%                 elems(m,:) = elems(m,idx);
%             end
%             
%             function idx = sortXYZ(vecs)
%                 hlp = vecs - min(vecs(:));
%                 dist = max(abs(vecs(:)))*100;
%                 hlp(3,:) = hlp(3,:)*dist^2;
%                 hlp(2,:) = hlp(2,:)*dist;
%                 [~, idx] = sort(sum(hlp,1),'ascend');
%             end
        end
    end
    
    methods(Access=private)
        
        function data = rawReadElems(~, filename)
            delimiter = ' ';

            %% Read columns of data as strings:
            % For more information, see the TEXTSCAN documentation.
            formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

            %% Open the text file.
            fileID = fopen(filename,'r');

            %% Read columns of data according to format string.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);

            %% Close the text file.
            fclose(fileID);

            %% Convert the contents of columns containing numeric strings to numbers.
            % Replace non-numeric strings with NaN.
            raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
            for col=1:length(dataArray)-1
                raw(1:length(dataArray{col}),col) = dataArray{col};
            end
            numericData = NaN(size(dataArray{1},1),size(dataArray,2));

            for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]
                % Converts strings in the input cell array to numbers. Replaced non-numeric
                % strings with NaN.
                rawData = dataArray{col};
                for row=1:size(rawData, 1);
                    % Create a regular expression to detect and remove non-numeric prefixes and
                    % suffixes.
                    regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
                    try
                        result = regexp(rawData{row}, regexstr, 'names');
                        numbers = result.numbers;

                        % Detected commas in non-thousand locations.
                        invalidThousandsSeparator = false;
                        if any(numbers==',');
                            thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                            if isempty(regexp(thousandsRegExp, ',', 'once'));
                                numbers = NaN;
                                invalidThousandsSeparator = true;
                            end
                        end
                        % Convert numeric strings to numbers.
                        if ~invalidThousandsSeparator;
                            numbers = textscan(strrep(numbers, ',', ''), '%f');
                            numericData(row, col) = numbers{1};
                            raw{row, col} = numbers{1};
                        end
                    catch me
                    end
                end
            end


            %% Replace non-numeric cells with NaN
            R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
            raw(R) = {NaN}; % Replace non-numeric cells

            %% Create output variable
            data = cell2mat(raw);
        end
        
        function [nodes, nodenrs] = rawReadNodes(~, filename)
            %% Start MatLab generated script
            delimiter = ' ';

            % Read columns of data as strings:
            formatSpec = '%*s%*s%*s%*s%s%*s%s%[^\n\r]';

            %% Open the text file.
            fileID = fopen(filename,'r');

            %% Read columns of data according to format string.
            % This call is based on the structure of the file used to generate this
            % code. If an error occurs for a different file, try regenerating the code
            % from the Import Tool.
            dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);

            %% Close the text file.
            fclose(fileID);

            %% Convert the contents of columns containing numeric strings to numbers.
            % Replace non-numeric strings with NaN.
            raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
            for col=1:length(dataArray)-1
                raw(1:length(dataArray{col}),col) = dataArray{col};
            end
            numericData = NaN(size(dataArray{1},1),size(dataArray,2));

            % Converts strings in the input cell array to numbers. Replaced non-numeric
            % strings with NaN.
            rawData = dataArray{2};
            for row=1:size(rawData, 1);
                % Create a regular expression to detect and remove non-numeric prefixes and
                % suffixes.
                regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
                try
                    result = regexp(rawData{row}, regexstr, 'names');
                    numbers = result.numbers;

                    % Detected commas in non-thousand locations.
                    invalidThousandsSeparator = false;
                    if any(numbers==',');
                        thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                        if isempty(regexp(thousandsRegExp, ',', 'once'));
                            numbers = NaN;
                            invalidThousandsSeparator = true;
                        end
                    end
                    % Convert numeric strings to numbers.
                    if ~invalidThousandsSeparator;
                        numbers = textscan(strrep(numbers, ',', ''), '%f');
                        numericData(row, 2) = numbers{1};
                        raw{row, 2} = numbers{1};
                    end
                catch me
                end
            end

            %% Split data into numeric and cell columns.
            rawNumericColumns = raw(:, 2);
            rawCellColumns = raw(:, 1);


            %% Replace non-numeric cells with NaN
            R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
            rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

            %% Allocate imported array to column variable names
            nodenrs = rawCellColumns(:, 1);
            nodes = cell2mat(rawNumericColumns(:, 1));
        end
    end
    
end

