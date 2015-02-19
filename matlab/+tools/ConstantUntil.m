classdef ConstantUntil < tools.AFunGen
    %CONSTANTUNTIL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access=private)
        stoptime;
        rampperc;
    end
    
    methods
        function this = ConstantUntil(time, percent)
            if nargin < 2
                percent = .005;
            end
            this.stoptime = time;
            this.rampperc = percent;
        end
        
        function [fhandle, df] = getFunction(this)
            st = this.stoptime;
            p = this.rampperc;
            rampstart = st*(1-p);
                          % 1 for t < rampstart  
            fhandle = @(t)(t < rampstart) + ...
                (t >= rampstart & t < st).*(1-(t-st*(1-p))/(st*p)); 
                % 1 to 0 over ramptime, zero after
            df = [];
        end
        
        function str = getConfigStr(this)
            str = sprintf('Endtime: %g, Ramp percent: %g',this.stoptime,this.rampperc);
        end
    end
    
end

