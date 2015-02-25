classdef QuadToLinear < tools.AFunGen
    % Returns the modified markert law functions for the OVERALL energy density
    % funcion derivative w.r.t. C (i.e. INCLUDING the 1/lam^2 prefactor
    % from chain rule!)
    
    properties(SetAccess=private)
        M;
        lam0;
    end
    
    methods
        function this = QuadToLinear(lam0, M)
            if nargin < 2
                M = 100;
                if nargin < 1
                    lam0 = 2;
                end
            end
            this.lam0 = lam0;
            this.M = M;
        end
        
        function [fhandle, dfhandle] = getFunction(this)
            m = this.M;
            l = this.lam0;
            a = m/2/(l-1);
            fl = a*l*l-2*a*l+a;
            fhandle = @(t)(t>=1 & t<l).*a.*(t.*t-2*t+1) + (t>=l).*(fl + (t-l)*m);
            dfhandle = @(t)(t>=1 & t<l).*a.*(2*t-2) + (t>=l)*m;
        end
        
        function str = getConfigStr(this)
            str = sprintf('lam0: %g, M: %g',this.lam0,this.M);
        end
        
        function plot(this, range, varargin)
            if nargin < 2
                range = [1 1.2*this.lam0];
            end
            plot@tools.AFunGen(this, range, varargin{:});
            if (range(2) > this.lam0)
                f = this.getFunction;
                ax = get(gcf,'Children');
                ax = ax(2); % second one is the left one - hope this is reproducible
                hold(ax,'on');
                plot(ax,this.lam0,f(this.lam0),'rx','MarkerSize',16);
            end
        end
    end
    
    methods(Static)
        function qfungen = fromLinearizedMarkert(mfungen)
            mfun = mfungen.getFunction;
            m = mfungen.max_modulus;
            t0 = 2*mfungen.t0;
            ft0 = mfun(t0);
            l0 = 2*(t0-ft0/m)-1;
            qfungen = tools.QuadToLinear(l0,m);
        end
    end
    
end

