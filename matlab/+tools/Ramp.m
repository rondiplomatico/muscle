classdef Ramp < tools.AFunGen
    
    properties(Access=private)
        ramptime;
        max;
        starttime;
    end
    
    methods
        function this = Ramp(ramptime, max, starttime)
            if nargin < 3
                starttime = 0;
                if nargin < 2
                    max = 1;
                    if nargin < 1
                        ramptime = 1;
                    end
                end
            end
            this.ramptime = ramptime;
            this.max = max;
            this.starttime = starttime;            
        end
        
        function fhandle = getFunction(this)
            rt = this.ramptime;
            if rt < 0
                fhandle = @(t)0;
            else
                maxval = this.max;
                start = this.starttime;
                fhandle = @(t)(t >= start) .* (maxval * (((t-start)<rt).*(t-start)/rt + (t>=rt+start)));
            end
        end
        
        function str = getConfigStr(this)
            str = sprintf('Start: %g, Ramptime: %g, Max:%g',this.starttime,this.ramptime,this.max);
        end
    end
    
end

