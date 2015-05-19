classdef CubicToLinear < tools.AFunGen
    % Returns the modified markert law functions for the OVERALL energy density
    % funcion derivative w.r.t. C (i.e. INCLUDING the 1/lam^2 prefactor
    % from chain rule!)
    
    properties(SetAccess=private)
        M;
        lam0;
    end
    
    methods
        function this = CubicToLinear(lam0, M)
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
            a = m/(l-1)^2; 
            fl = a*(l.^3/3-l.^2+l-1/3);
            fhandle = @(t)(t>=1 & t<l).*a.*(t.^3/3-t.^2+t-1/3) + (t>=l).*(fl + (t-l)*m);
            dfhandle = @(t)(t>=1 & t<l).*a.*(t.^2-2*t+1) + (t>=l)*m;
        end
        
        function str = getConfigStr(this)
            str = sprintf('lam0: %g, M: %g',this.lam0,this.M);
        end
        
        function pm = plot(this, varargin)
            % @fixme
            i = inputParser;
            i.KeepUnmatched = true;
            i.addParamValue('R',[1 1.2*this.lam0]);
            i.parse(varargin{:});
            res = i.Results;
            pm = plot@tools.AFunGen(this, varargin{:});
            if (res.R(2) > this.lam0)
                mc = metaclass(this);
                for k = 1:length(pm.Figures)
                    if ishandle(pm.Figures(k)) && strcmpi(get(pm.Figures(k),'Tag'),...
                            strrep(mc.Name,'.','_'))
                        f = this.getFunction;
                        ax = get(pm.Figures(k),'Children');
                        hold(ax(1),'on');
                        plot(ax(1),this.lam0,f(this.lam0),'rx','MarkerSize',16);
                        break;
                    end
                end
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
            qfungen = tools.CubicToLinear(l0,m);
        end
    end
    
end

