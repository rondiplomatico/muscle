classdef Sinus < tools.AFunGen
    
    properties(Access=private)
        freq;
        foffset;
        voffset;
    end
    
    methods
        function this = Sinus(freq, foffset, voffset)
            if nargin < 3
                voffset = 0;
                if nargin < 2
                    foffset = 0;
                    if nargin < 1
                        freq = 3;
                    end
                end
            end
            this.freq = freq;
            this.foffset = foffset;
            this.voffset = voffset;
        end
        
        function fhandle = getFunction(this)
            f = this.freq;
            fo = this.foffset;
            vo = this.voffset;
            fhandle = @(t)sin(t/1000*f*2*pi+fo)+vo;
        end
        
        function str = getConfigStr(this)
            str = sprintf('Frequency: %g [Hz], ArgOffset: %g, valueOffset: %g',...
                this.freq,this.foffset,this.voffset);
        end
    end
    
end

