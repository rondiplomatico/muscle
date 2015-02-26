classdef MarkertLawOriginal < tools.AFunGen
    % Returns the original markert law functions for the OVERALL energy density
    % funcion derivative w.r.t. C (i.e. INCLUDING the 1/lam^2 prefactor
    % from chain rule!)
    
    properties(SetAccess=private)
        b;
        d;
    end
    
    methods
        function this = MarkertLawOriginal(b, d)
            if nargin < 2
                d = 16.6;
                if nargin < 1
                    b = 7990;
                end
            end
            this.b = b;
            this.d = d;
        end
        
        function [fhandle, dfhandle, fbdhandle, dfbdhandle] = getFunction(this)
            b1 = this.b;
            d1 = this.d;            
            
            %% Regular function
%             fhandle = @(t)b1./t.^2.*(t.^d1-1);
%             dfhandle = @(t)b1./t.^3.*((d1-2).*t.^d1+2);
%             fbdhandle = @(t,b1,d1)b1./t.^2.*(t.^d1-1);
%             dfbdhandle = @(t,b1,d1)b1./t.^3.*((d1-2).*t.^d1+2);
            
            %% Variant with support [1, \infty]
            fhandle = @(t)max(0,b1./t.^2.*(t.^d1-1));
            dfhandle = @(t)(t>1).*(b1./t.^3.*((d1-2).*t.^d1+2));
            fbdhandle = @(t,b1,d1)max(0,b1./t.^2.*(t.^d1-1));
            dfbdhandle = @(t,b1,d1)max(0,b1./t.^3.*((d1-2).*t.^d1+2));
            
            %% Variant without dlam/dC term (1/2\lambda)
%             fhandle = @(t)b1./t.*(t.^d1-1);
%             dfhandle = @(t)b1./t.^2.*((d1-2).*t.^d1+2);
%             fbdhandle = @(t,b1,d1)b1./t.*(t.^d1-1);
%             dfbdhandle = @(t,b1,d1)b1./t.^2.*((d1-2).*t.^d1+2);
        end
        
        function str = getConfigStr(this)
            str = sprintf('b: %g, d: %g',this.b,this.d);
        end
        
        function pm = plot(this, range, varargin)
            if nargin < 2
                range = [.9 1.2];
            end
            pm = plot@tools.AFunGen(this, range, varargin{:});
        end
    end
    
    methods(Static)
        function createScratchPlot
            range = [.9 1.4];
            so = tools.MarkertLawOriginal(7990,16.6);
            pm = so.plot(range,'b');
            ax = get(gcf,'Children');
            ax = flipud(ax);
            hold(ax(1),'on');
            hold(ax(2),'on');
            ylim(ax(1),[0 16500]);
            ylim(ax(2),[0 4.1e5]);
            s = tools.MarkertLaw(7990,16.6,1e9);
            s.plot(range,ax,'r');
            s = tools.MarkertLaw(1000*7990,4.5,1e9);
            s.plot(range,ax,'g');
            legend(ax(1),'Original b=7990,d=16.6','Modified b=7990,d=16.6','Modified b=7990000,d=4.5');%,'Location','NorthWest');
            pm.savePlots('/home/dwirtz/software/muscle/latex/img','Format','eps');
        end
    end
    
end

