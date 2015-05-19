classdef AFunGen < handle
    %AFUNGEN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xLabel = 't [ms]';
        yLabel = 'value';
    end
    
    methods
        function pm = plot(this, varargin)
            [f, df] = this.getFunction;
            args = {};
            if ~isempty(df)
                args = {false, 1, 2};
            end
            pm = PlotManager(args{:});
            pm.UseFileTypeFolders = false;
            pm.LeaveOpen = true;
            
            i = inputParser;
            i.KeepUnmatched = true;
            i.addParamValue('PM',pm,@(v)isa(v,'PlotManager'));
            i.addParamValue('R',0:.1:1000);
            i.addParamValue('P',{});
            i.addParamValue('AX',[],@(v)all(ishandle(v)));
            i.parse(varargin{:});
            res = i.Results;
            pm = res.PM;
            if length(res.R) == 2
                range = linspace(res.R(1),res.R(2),2000);
            else
                range = res.R;
            end
            
            if ~isempty(res.AX)
                hold(res.AX(1),'on');
                plot(res.AX(1),range,f(range),varargin{:});
                hold(res.AX(2),'on');
                plot(res.AX(2),range,f(range),varargin{:});
            else
                mc = metaclass(this);
                ax = pm.nextPlot(mc.Name,...
                    sprintf('Plot of %s\n%s',mc.Name,this.getConfigStr),...
                    this.xLabel,this.yLabel);
                plot(ax,range,f(range),res.P{:});
                if ~isempty(df)
                    ax = pm.nextPlot([mc.Name '_deriv'],...
                        sprintf('Plot of %s-derivative\n%s',mc.Name,this.getConfigStr),...
                    this.xLabel,[this.yLabel '/' this.xLabel]);
                    plot(ax,range,df(range),res.P{:});
                end
            end
            if nargout > 0
                pm.done;
            end
        end
        
        function plottofigure(this, fignr, range, varargin)
            if nargin < 3
                range = 0:.1:1000;
                if nargin < 2
                    fignr = gcf;
                end
            elseif length(range) == 2
                range = linspace(range(1),range(2),2000);
            end
            ax = get(fignr,'Children');
            if numel(ax) < 2
                error('Need two axes handles.');
            end
            this.plot(range, ax, varargin{:});
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

