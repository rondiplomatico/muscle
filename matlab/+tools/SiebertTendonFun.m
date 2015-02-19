classdef SiebertTendonFun < tools.AFunGen
    % Returns the modified markert law functions for the OVERALL energy density
    % funcion derivative w.r.t. C (i.e. INCLUDING the 1/lam^2 prefactor
    % from chain rule!)
    
    properties(SetAccess=private)
        l0;
        lm = .742; % [mm]
        lam_m;
        F1 = 1345; % [mN]
    end
    
    methods
        function this = SiebertTendonFun(l0)
            if nargin < 1
                l0 = [];
            end
            if ~isempty(l0)
                this.l0 = l0;
                this.lam_m = 1+this.lm/l0;
            end
        end
        
        function [fhandle, dfhandle] = getFunction(this)
            k = 4.36;
            F_1 = this.F1;
            c = F_1/(exp(k)-1);
            K = 8030; % [mN/mm]
            
            lnull = this.l0;
            if ~isempty(lnull)
                lm = this.lam_m;
                fhandle = @(l)(l<=lm)*c.*(exp(k*(l-1)/(lm-1))-1) + (l > lm).*(F_1+K*(l-lm)*lnull);
            else
                lm = .742;
                fhandle = @(l)(l<=lm)*c.*(exp(k*l/lm)-1) + (l > lm).*(F_1+K*(l-lm));
            end
            dfhandle = [];
        end
        
        function str = getConfigStr(this)
            if ~isempty(this.l0)
                str = sprintf('l0: %g, \\lambda_m: %g',this.l0,this.lam_m);
            else
                str = sprintf('Original from paper Siebert 2012');
            end
        end
        
        function plot(this, range)
            point = this.lm;
            if ~isempty(this.l0)
                point = this.lam_m;
                if nargin < 2
                    range = [1 1.2*point];
                end
                lab = '\lambda [-]';
            else
                if nargin < 2
                    range = [0 4];
                end
                lab = '\Delta stretch [mm]';
            end
            plot@tools.AFunGen(this, range);
            
%             ax = get(gcf,'Children');
%             ax = ax(2); % second one is the left one - hope this is reproducible
            
            ax = gca;
            hold(ax,'on');
            plot(ax,point,this.F1,'rx','MarkerSize',16);
            xlabel(ax,lab);
            ylabel(ax,'force [mN]');
        end
    end
    
end

