classdef AFunGen < handle
    %AFUNGEN Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function plot(this, range)
            if nargin < 2
                range = 0:.1:1000;
            elseif length(range) == 2
                range = linspace(range(1),range(2),2000);
            end
            [f, df] = this.getFunction;
            args = {};
            if ~isempty(df)
                args = {false, 1, 2};
            end
            pm = PlotManager(args{:});
            pm.LeaveOpen = true;
            mc = metaclass(this);
            ax = pm.nextPlot(mc.Name,sprintf('Plot of %s\n%s',mc.Name,this.getConfigStr),'t [ms]','value');
            plot(ax,range,f(range));
            if ~isempty(df)
                ax = pm.nextPlot(mc.Name,sprintf('Plot of %s-derivative\n%s',mc.Name,this.getConfigStr),'t [ms]','value');
                plot(ax,range,df(range));
            end
            pm.done;
        end
        
        function sum = plus(this, other)
            % Provides an override for the simple sum of two AFunGen
            if ~isa(this,'tools.AFunGen') || ~isa(other,'tools.AFunGen')
                error('Addition not defined for non-AFunGen classes.');
            end
            sum = tools.FuncSum(this, other);
        end
    end
    
    methods(Abstract)
        [fhandle, dfhandle] = getFunction(this);
        
        str = getConfigStr(this);
    end
    
end

