classdef MarkertLaw < tools.AFunGen
    % Returns the modified markert law functions for the OVERALL energy density
    % funcion derivative w.r.t. C (i.e. INCLUDING the 1/lam^2 prefactor
    % from chain rule!)
    
    properties(Constant)
        LatestLinearizationLambda = 5;
    end
    
    properties(SetAccess=private)
        b;
        d;
        max_modulus;
        t0 = [];
        ft0 = [];
    end
    
    methods
        function this = MarkertLaw(b, d, max_modulus)
            if nargin < 3
                max_modulus = [];
                if nargin < 2
                    d = 2;
                    if nargin < 1
                        b = 1;
                    end
                end
            end
            this.b = b;
            this.d = d;
            if ~isempty(max_modulus) && d > 1
                [f,df] = this.getFunction;
                rightval = 2;
                found = false;
                while ~found
                    pts = linspace(1,rightval,3);
                    % Compute maxit to the numer of iterations where the
                    % interval length divided by 2^maxit < eps
                    maxit = ceil(log((pts(3)-pts(1))/eps)/log(2));
                    it = 0; diff = Inf;
                    while diff > sqrt(eps) && it < maxit;
                        fpts = df(pts)-max_modulus;
                        if fpts(2) < 0
                            pts = [pts(2) (pts(2)+pts(3))/2 pts(3)];
                        else
                            pts = [pts(1) (pts(1)+pts(2))/2 pts(2)];
                        end
                        diff = abs(fpts(2));
                        it = it+1;
                        if pts(2) > this.LatestLinearizationLambda
                            found = true;
                            break;
                        end
                    end
                    if rightval-pts(2) < 1e-10
                        rightval = rightval*2;
                        %fprintf('Increasing right interval value to %d\n',rightval);
                    else
                        found = true;
                    end
                end
                this.t0 = pts(2);
                fprintf('Found t0=%5.3g with diff %g to max modulus/derivative of %g after %d iterations (b=%g, d=%g).\n',...
                    this.t0,diff,max_modulus,it,b,d);
                this.ft0 = f(this.t0);
                this.max_modulus = max_modulus;
            end
        end
        
        function [fhandle, dfhandle, fbdhandle, dfbdhandle] = getFunction(this)
            b1 = this.b;
            d1 = this.d;
            tlin = this.t0;
            if ~isempty(tlin)
                mm = this.max_modulus;
                ftlin = this.ft0;
                %fhandle = @(t)(t<tlin).*max(0,b1.*(t.^d1-1)) + (t>=tlin).*(ftlin + mm*(t-tlin));
                %dfhandle = @(t)(t<tlin).*max(0,b1.*((d1-2).*t.^d1+2)) + (t>=tlin)*mm;
                fhandle = @(t)(t >= 1 & t<tlin).*(b1./t.^2.*(t.^d1-1).*(t-1).^2)...
                    + (t>=tlin).*(ftlin + mm*(t-tlin));
                dfhandle = @(t)(t >= 1 & t<tlin).*(b1.*(t-1)./t.^3.*((d1*(t-1)+2).*t.^d1-2)) ...
                    + (t>=tlin)*mm;
                fbdhandle = @(t,b1,d1)(t >= 1 & t<tlin).*(b1./t.^2.*(t.^d1-1).*(t-1).^2)...
                    + (t>=tlin).*(ftlin + mm*(t-tlin));
                dfbdhandle = @(t,b1,d1)(t >= 1 & t<tlin).*(b1.*(t-1)./t.^3.*((d1*(t-1)+2).*t.^d1-2)) ...
                    + (t>=tlin)*mm;
            else
%                 fhandle = @(t)max(0,b1.*(t.^d1-1));
%                 dfhandle = @(t)max(0,b1.*((d1-2).*t.^d1+2));
                fhandle = @(t)(t>=1).*(b1./t.^2.*(t.^d1-1).*(t-1).^2);
                dfhandle = @(t)(t>=1).*b1.*(t-1)./t.^3.*((d1*(t-1)+2).*t.^d1-2);
                fbdhandle = @(t,b1,d1)(t>=1).*(b1./t.^2.*(t.^d1-1).*(t-1).^2);
                dfbdhandle = @(t,b1,d1)(t>=1).*b1.*(t-1)./t.^3.*((d1*(t-1)+2).*t.^d1-2);
            end
        end
        
        function str = getConfigStr(this)
            str = sprintf('b: %g, d: %g, MaxModulus:%g',this.b,this.d,this.max_modulus);
        end
        
        function plot(this, range, varargin)
            if nargin < 2
                range = [1 1.2*this.t0];
            end
            plot@tools.AFunGen(this, range, varargin{:});
            if (range(2) > this.t0)
                ax = get(gcf,'Children');
                ax = ax(2); % second one is the left one - hope this is reproducible
                hold(ax,'on');
                plot(ax,this.t0,this.ft0,'rx','MarkerSize',16);
            end
        end
    end
    
    methods(Static)
        function test_Linarization
            C = Utils.createCombinations([1 5],linspace(3,40,5),linspace(100,1e6,5));
            x = linspace(.3,6,300);
            pm = PlotManager(false,1,2);
            pm.LeaveOpen = true;
            ax = pm.nextPlot('markertfun_lintest','Plots for various parameters','lambda','value');
            axl = pm.nextPlot('markertfun_lintest_log','Plots for various parameters','lambda','value');
            for i = 1:size(C,2)
                ml = tools.MarkertLaw(C(1,i),C(2,i),C(3,i));
%                 ml.plot(x);                
                [f,df] = ml.getFunction;
                fx = f(x);
                plot(ax,x,fx,'b');
                semilogy(axl,x,fx,'b');
                fprintf('t0=%g, f(t0)=%g for %s\n',ml.t0, ml.ft0, ml.getConfigStr);
                
                hold(ax,'on');
                hold(axl,'on');
                
                ml = tools.MarkertLaw(C(1,i),C(2,i),[]);
%                 ml.plot(x);
                [f,df] = ml.getFunction;
                fx = f(x);
                plot(ax,x,fx,'r');
                semilogy(axl,x,fx,'r');
            end
            pm.done;
            ylim(ax,[0 1e6]);
        end
    end
    
end

