classdef AFunGen < handle
    %AFUNGEN Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function plot(this, range)
            if nargin < 2
                range = 0:.1:1000;
            elseif length(range) == 2
                range = linspace(range(1),range(2),1000);
            end
            f = this.getFunction;
            pm = PlotManager;
            pm.LeaveOpen = true;
            mc = metaclass(this);
            ax = pm.nextPlot(mc.Name,sprintf('Plot of %s\n%s',mc.Name,this.getConfigStr),'t [ms]','value');
            plot(ax,range,f(range));
            pm.done;
        end
    end
    
    methods(Abstract)
        fhandle = getFunction(this);
        
        str = getConfigStr(this);
    end
    
end

