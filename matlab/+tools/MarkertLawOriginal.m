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
                d = 2;
                if nargin < 1
                    b = 1;
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
        
        function pm = plot(this, range)
            if nargin < 2
                range = [.5 2];
            end
            pm = plot@tools.AFunGen(this, range);
        end
    end
    
end

