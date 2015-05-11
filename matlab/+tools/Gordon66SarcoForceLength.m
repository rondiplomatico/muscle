classdef Gordon66SarcoForceLength < tools.PiecewiseLinear
    
    properties(Access=private)
        l0;
        cdata;
    end
    
    methods
        function this = Gordon66SarcoForceLength(l0)
            if nargin < 1
                l0 = 2.05;
            end
            
            % Data from Paper Gordon1966
            dopt_tab = [1.98 97.09; 2.03 98.6; 2.08 99.87; 2.13 100; 2.18 100.09; 2.23 98.85; 2.28 97.06];
            % Roughly "estimated"
            dup = [1.3 7; 1.4 25; 1.45 35; 1.5 49; 1.6 68; 1.65 78; 1.7 85; 1.75 88; 1.8 91; 1.9 95];
            ddow = [2.3 96.6; 2.8 60; 3.0 45; 3.4 20; 3.6 6];
            curve = [dup; dopt_tab; ddow];
            
            % Sort
            %[~, idx] = sort(curve(:,1));
            %curve = curve(:,idx);
            
            % Normalize
            curve(:,2) = curve(:,2)/max(curve(:,2));
            
            % Pad zeros
            curve = [1.1 0; 1.25 0; curve; 3.7 0; 3.8 0];
            
            sel = [2 9 13 19 25];
            this = this@tools.PiecewiseLinear(curve(sel,1)',curve(sel,2)');
            this.l0 = l0;
            this.cdata = curve;
        end
        
        function str = getConfigStr(this)
            str = sprintf('L_0=%g',this.l0);
        end
        
        function plot(this, varargin)
            pm = plot@tools.PiecewiseLinear(this, varargin{:});
            
            ax = pm.nextPlot('gordon66_fit','Splint fits','length','force');
            hold(ax,'on');
            
            [len, proc] = this.transform(this.cdata(:,1)',this.cdata(:,2)');
            
            plot(len,proc,'r-x');
            
            % PChip/Spline plots
            slen = linspace(min(len)*.99,max(len)*1.01,200);
            sproc = spline(len,proc,slen);
            pp = pchip(len,proc);
            pproc = ppval(slen,pp);
            plot(ax,slen,sproc,'k-',slen,pproc,'m-');
            
            pol = polyfit(len,proc,8);
            %fl = @(x)pol(1)*x.*x.*x + pol(2)*x.*x + pol(3)*x + pol(4);
            plot(ax,slen,polyval(pol,slen),'g');
            
            axis(ax,'tight');
        end
        
    end
    
    methods(Access=protected)
        function [x,y] = transform(this, x, y)
            x  = x / this.l0;
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