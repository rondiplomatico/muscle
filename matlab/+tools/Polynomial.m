classdef Polynomial < tools.AFunGen
    
    properties(Access=private)
        poly;
    end
    
    methods
        function this = Polynomial(poly)
            if nargin < 1
                poly = 1;
            end
            this.poly = poly;
        end
        
        function [fhandle, dfhandle] = getFunction(this)
            p = this.poly;
            fhandle = @(t)polyval(p,t);
            dp = (8:-1:1) .* p(1:end-1);
            dfhandle = @(t)polyval(dp,t);
        end
        
        function str = getConfigStr(this)
            str = sprintf('Degree=%d',length(this.poly));
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