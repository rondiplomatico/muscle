classdef Gordon66SarcoForceLength < tools.AFunGen
    
    properties(Access=private)
        l0;
        cdata;
    end
    
    methods
        function this = Gordon66SarcoForceLength(l0)
            if nargin < 1
                l0 = 2.05;
            end
            this.l0 = l0;
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
            
            this.cdata = curve;
        end
        
        function [fhandle, dfhandle] = getFunction(this)
            curve = this.cdata;
            % Scale to l_0 (standard/resting length)
            curve(:,1) = curve(:,1) / this.l0;
            
            % Extract piecewise linear function
            parts = [2 9 13 19 25];
            np = length(parts)-1;
            funs = cell(1,np);
            dfuns = cell(1,np);
            for k=1:np
                xy1 = curve(parts(k),:);
                xy2 = curve(parts(k+1),:);
                cond = sprintf('(t>=%g & t <%g)',xy1(1),xy2(1));
                a = (xy2(2)-xy1(2))/(xy2(1)-xy1(1));
                b = xy1(2)-a*xy1(1);
                fun = sprintf('(%g*t+%g)',a,b);
                funs{k} = [cond '.*' fun];
                dfuns{k} = [cond '.*' sprintf('%g',a)];
            end
            funstr = sprintf('@(t)%s',Utils.implode(funs,' + '));
            dfunstr = sprintf('@(t)%s',Utils.implode(dfuns,' + '));
            fhandle = eval(funstr);
            dfhandle = eval(dfunstr);
        end
        
        function str = getConfigStr(this)
            str = sprintf('L_0=%g',this.l0);
        end
        
        function plot(this, range)
            curve = this.cdata;
            % ALSO Scale to l_0 (standard/resting length)
            curve(:,1) = curve(:,1) / this.l0;
            
            if nargin < 2
                range = linspace(min(curve(:,1))*0.95,max(curve(:,1))*1.05,1000);
            end
            plot@tools.AFunGen(this, range);
            
            ax = gca;
            hold(ax,'on');
            
            len = curve(:,1)';
            proc = curve(:,2)';
            
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