classdef PiecewiseLinear < tools.AFunGen
    
    properties(Access=protected)
        pts;
        vals;
    end
    
    methods
        function this = PiecewiseLinear(x,y)
            if nargin < 2
                y = x;
                x = 1:size(x,2);
                if nargin < 1
                    x = 0:2;
                    y = [0 1 0];
                end
            end
            if length(x) ~= length(unique(x))
                error('Have to provide disjoint x locations');
            end
            [x, idx] = sort(x,'ascend');
            y = y(idx);
            this.pts = x;
            this.vals = y;
        end
        
        function [fhandle, dfhandle] = getFunction(this)
            [x,y] = this.transform(this.pts, this.vals);
            xy = [x;y];
            nf = size(xy,2)-1;
            funs = cell(1,nf);
            dfuns = cell(1,nf);
            for k=1:nf
                xy1 = xy(:,k);
                xy2 = xy(:,k+1);
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
            str = sprintf('Points: %d',size(this.pts,2));
        end
        
        function pm = plot(this, range)
            if nargin < 2
                x = this.transform(this.pts,this.vals);
                range = [min(x) max(x)];
            end
            pm = plot@tools.AFunGen(this, range);
        end
        
    end
    
    methods(Access=protected)
        function [x,y] = transform(this, x, y)
           % gives the chance for subclasses to manipulate the x,y data
           % before the piecewise nonlinear function is created
           %
           % do nothing by default
        end
    end
    
    methods(Static)
        function test_PiecewiseLinear
            p = tools.PiecewiseLinear;
            p.plot;
            p = tools.PiecewiseLinear(rand(1,7));
            p.plot;
            p = tools.PiecewiseLinear(rand(1,7),rand(1,7));
            p.plot;
        end
    end
    
end