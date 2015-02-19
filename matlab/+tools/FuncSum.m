classdef FuncSum < tools.AFunGen
    
    properties(Access=private)
        fungens;
    end
    
    methods
        function this = FuncSum(varargin)
            if length(varargin) < 2
                error('Must have at least two summands!');
            end
            this.fungens = varargin;
        end
        
        function [fhandle, dfhandle] = getFunction(this)%#ok
            f = this.fungens;
            names = cell(1,length(f));
            for k = 1:length(f)
                mc = metaclass(f{k});
                [~,name] = fileparts(which(mc.Name));
                names{k} = sprintf('%s%d',name,k);
                eval(sprintf('%s = f{%d}.getFunction;',names{k},k));
            end
            % Add (t) argument
            names = cellfun(@(s)sprintf('%s(t)',s),names,'UniformOutput',false);
            % Join with +
            sumstr = Utils.implode(names,' + ');
            % Eval to readable handle
            eval(sprintf('fhandle = @(t)%s;',sumstr));
            dfhandle = [];
        end
        
        function str = getConfigStr(this)
            str = this.fungens{1}.getConfigStr;
            for k = 2:length(this.fungens)
                str = sprintf('/ %s',this.fungens{k}.getConfigStr);
            end
        end
    end
    
end

% Speed test
% all(1:3) = 0;
% for k=1:1000
%     slen = min(len)*.99 + rand*(max(len)*1.01-min(len)*.99);
%     t = tic;
%     %sproc = spline(len,proc,slen);
%     sproc = pchip(len,proc,slen);
%     all(1) = all(1) + toc(t);
%     l = rand*3;
%     t = tic;
%     fl(l);
%     all(2) = all(2) + toc(t);
%     t = tic;
%     polyval(pol,slen);
%     all(3) = all(3) + toc(t);
% end
% all / 1000