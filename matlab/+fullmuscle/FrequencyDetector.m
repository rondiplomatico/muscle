classdef FrequencyDetector < KerMorObject
    % FrequencyDetector:
    % 
    % The old/original version was made by @ts
    %
    % @author Daniel Wirtz @date 2014-09-23
    %
    % @new{0,7,dw,2014-09-23} Added this class.
    
    properties
        PeakOnThreshold = 40;
        PeakOffThreshold = 15;
    end
    
    properties(SetAccess=private)
        % The current frequency
        Frequency = 0;
        
        % The number of concurrently tracked peaks
        WindowSize = 4;
    end
    
    properties(Access=private)
        % times of last p peaks
        peaktimes;
        
        % flag to remember if currently in peak
        ispeak;
        
        % The number of concurrently processed signals
        ns;
    end
    
    methods
        function this = FrequencyDetector(num_signals,windowsize,sig_on,sig_off)
            if nargin < 1
                num_signals = 1;
            end
            this.ns = num_signals;
            if nargin > 1
                this.WindowSize = windowsize;
                if nargin > 2
                    this.PeakOnThreshold = sig_on;
                    if nargin > 3
                        this.PeakOffThreshold = sig_off;
                    end
                end
            end
            this.reset;
        end
        
        function reset(this)
            this.Frequency = zeros(1,this.ns);
            this.peaktimes = -Inf(this.WindowSize,this.ns);
            this.ispeak = false(1,this.ns);
        end
        
        function s = processSignal(this, t, sig)
            % A new peak is detected
            peak_on = ~this.ispeak & (sig > this.PeakOnThreshold);
            peak_off = this.ispeak & sig < this.PeakOffThreshold;
            if any(peak_on)
                % We're on peak!
                this.ispeak(peak_on) = true;
                % Get current frequency for windowsize peaks over time
                this.Frequency(peak_on) = this.WindowSize ./ (t-this.peaktimes(1,peak_on));
                % "Remove" last peak and Shift peaktimes on
                this.peaktimes(1:end-1,peak_on) = this.peaktimes(2:end,peak_on);
                % Save current peak time
                this.peaktimes(end,peak_on) = t;
            end
            this.ispeak(peak_off) = false;
            s = 0;
        end
    end
    
    methods(Static)
        function test_FrequencyDetector
            fd = fullmuscle.FrequencyDetector(1,4);
            fdw3 = fullmuscle.FrequencyDetector(1,3);
            
            % Plain run
            m = models.motoneuron.Model;
            m.T = 500;
            [t,y] = m.simulate;
            
            freq = zeros(size(t));
            freqw3 = freq;
            for k=1:length(t)%#ok
                fd.processSignal(t(k),y(2,k));
                freq(k) = fd.Frequency;
                fdw3.processSignal(t(k),y(2,k));
                freqw3(k) = fdw3.Frequency;
            end
            plot(t,freq,'r',t,freqw3,'b');
            legend('Window=4','Window=3');
            
            % Integrated ODE-run
            freq_ode = zeros(size(t));
            fd.reset;
            k=1;
            slv = m.ODESolver;
            slv.odeopts = odeset(slv.odeopts,...
                'OutputFcn',@fd_update,...
                'OutputSel',2);
            m.simulate;
            plot(t,freq_ode);
            
            function s = fd_update(t,y,~)
                fd.processSignal(t,y);
                freq_ode(k) = fd.Frequency;
                k = k+1;
                s = 0;
            end
            
            % Multi-signal run
            fd = fullmuscle.FrequencyDetector(3);
            slv.odeopts = odeset(slv.odeopts,...
                'OutputFcn',[],...
                'OutputSel',[]);
            [~,y2] = m.simulate([0; 3]);
            [~,y3] = m.simulate([1; 7]);
            freq = zeros(3,length(t));
            for j=1:length(t)
                fd.processSignal(t(j),[y(2,j) y2(2,j) y3(2,j)]);
                freq(:,j) = fd.Frequency;
            end
            plot(t,freq');
        end
    end
    
end
