classdef MusclePlotter < muscle.MusclePlotter
% MusclePlotter: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2014-09-25
%
% @new{0,7,dw,2014-09-25} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    methods
        
        function this = MusclePlotter(sys)
            this = this@muscle.MusclePlotter(sys);
        end
        
        function [pm, h_geo] = plot(this, t, y, varargin)
            opts = this.parsePlotArgs(varargin);
            
            mc = this.Config;
            sys = this.System;
            sel = opts.MU;
            
            if isempty(opts.PM)
                if opts.GeoOnly
                    pm = PlotManager;
                else
                    pm = PlotManager(false,2,3);
                end
                pm.LeaveOpen = true;
            else
                pm = opts.PM;
            end
            
            [pd, t, y] = this.updatePlotData(this.plotdata, opts, t, y);
            this.plotdata = pd;
            
            if opts.Vid
                vw = VideoWriter(opts.Vid);
                vw.FrameRate = 25;
                vw.open;
            end
            
            if opts.Geo
    %             pmgeo = PlotManager;
    %             pmgeo.LeaveOpen = true;
                pmgeo = pm;
                h_geo = pmgeo.nextPlot('geo',sprintf('Deformation at t=%g',t(end)),'x [mm]','y [mm]');
                zlabel(h_geo,'z [mm]');
                axis(h_geo, this.getPlotBox(y));
                daspect([1 1 1]);
                view(h_geo, this.GeoView);
                hold(h_geo,'on');
            end
            
            if ~opts.GeoOnly
                
                if opts.Moto
                    h1 = pm.nextPlot('signal','Motoneuron signal','t [ms]','V_m');
                    axis(h1,[0 t(end) min(pd.moto_vm(:)) max(pd.moto_vm(:))]);
                    hold(h1,'on');
                end

                if opts.Sarco
%                     h2 = pm.nextPlot('force','Action potential','t [ms]','V_m');
%                     axis(h2,[0 t(end) min(pd.sarco_pot(:)) max(pd.sarco_pot(:))]);
%                     hold(h2,'on');

                    h3 = pm.nextPlot('force','Motorunit force signals','t [ms]','A_s');
                    val = pd.sarco_force;
                    axis(h3,[0 t(end) min(val(:)) max(val(:))]);
                    hold(h3,'on');
                    
                    h3b = pm.nextPlot('force','Weighted Activation at Elem1,GP1','t [ms]','A_s');
                    val = mc.FibreTypeWeights(1,:,1)*pd.sarco_force;
                    axis(h3b,[0 t(end) min(val(:)) max(val(:))]);
                    hold(h3b,'on');
                end

                if opts.Freq
                    hfreq = pm.nextPlot('frequency','Motoneuron frequency','t [ms]','Frequency [Hz]');
                    m = min(pd.freq(:));
                    M = max(pd.freq(:));
                    if opts.FreqDet && ~sys.f.UseFrequencyDetector
                        m = min(m, min(pd.freq_det(:)));
                        M = max(M, max(pd.freq_det(:)));
                    end
                    axis(hfreq,[0 t(end) m M]);
                    hold(hfreq,'on');
                end

                spos = mc.SpindlePositions;
                if opts.Spin
                    h_spin_l = pm.nextPlot('spindle_lambda',...
                        'Spindle lambda stretch','t [ms]','Lambda [L_0 = 1]');
                    axis(h_spin_l,[0 t(end) min(pd.spindle_lambda(:)) max(pd.spindle_lambda(:))]);
                    hold(h_spin_l,'on');

                    if opts.Aff
                        h4 = pm.nextPlot('spindle','Afferents','t [ms]','aff');
                        axis(h4,[0 t(end) min(pd.afferents(:)) max(pd.afferents(:))]);
                        hold(h4,'on');
                        affsel = [2*(sel-1)+1; 2*(sel-1)+2];        
                    end
                end
                
                h5 = pm.nextPlot('mean_current','Motoneuron input mean current','t [ms]','mean current');
                axis(h5,[0 t(end) min(pd.eff_mean_current(:)) max(pd.eff_mean_current(:))+eps]);
                hold(h5,'on');

                if opts.Ext
                    hext = pm.nextPlot('ext_neumann',...
                        'External pressure','t [ms]','Normal pressure [kPa]');
                    axis(hext,[0 t(end) min(pd.uneum(:)) max(pd.uneum(:))+eps]);
                    hold(hext,'on');

    %                 val = pd.uneum * pd.forcefactor;
    %                 hext_f = pm.nextPlot('ext_neumann',...
    %                     'External force','t [ms]','Normal force [mN]');
    %                 axis(hext_f,[0 t(end) min(val(:)) max(val(:))]);
    %                 hold(hext_f,'on');
                end
                
            end
                
            pm.done;
            
            fh = gcf;
            for ts = 1:length(t)
                % Quit if figure has been closed
                if ~ishandle(fh)
                    break;
                end
                
                if opts.Geo
                    yf = pd.yfull(:,ts);
                    this.plotGeometry(h_geo, t(ts), yf, ts, opts);
                    
                    if ~opts.GeoOnly && opts.Spin && ~isempty(spos)
                        % Plot spindle locations
                        for k = sel
                            u = yf(sys.idx_u_glob_elems(:,:,spos(1,k)));
                            spindle_pt = u*pd.Ngp(:,spos(2,k));
                            plot3(h_geo,spindle_pt(1),spindle_pt(2),...
                                spindle_pt(3),'.',...
                                'MarkerSize',20,'Color',[.3 1 .3]);
                        end 
                    end
                end
                
                if ~opts.GeoOnly
                
                    time_part = t(1:ts);
                    if opts.Moto
                        cla(h1);
                        plot(h1,time_part,pd.moto_vm(sel,1:ts));
                    end

                    if opts.Sarco
%                         cla(h2);
%                         plot(h2,time_part,pd.sarco_pot(sel,1:ts));

                        cla(h3);
                        force = pd.sarco_force(sel,1:ts);
                        plot(h3,time_part,force,'r');
                        
                        cla(h3b);
                        walpha = mc.FibreTypeWeights(1,sel,1) * force;
                        plot(h3b,time_part,walpha,'b');
                    end

                    cla(h5);
                    plot(h5,time_part,pd.ext_mean_current(1:ts),'b--');
                    plot(h5,time_part,pd.eff_mean_current(sel,1:ts),'r');
                    
                    if opts.Spin
                        cla(h_spin_l);
                        plot(h_spin_l, time_part, pd.spindle_lambda(sel,1:ts));
                        
                        % Also add the spindle mean current
                        plot(h5,time_part,pd.spindle_mean_current(1:ts),'g--');

                        if opts.Aff
                            cla(h4);
                            plot(h4,time_part,pd.afferents(affsel(:),1:ts)');
                            

    %                         cla(h6);
    %                         plot(h6,time_part,pd.spindle_single_mean_current(sel,1:ts));
    %                         plot(h6,time_part,pd.eff_mean_current(sel,1:ts),'r--');
                        end
                    end
%                     axis(h5,'tight');

                    if opts.Freq
                        cla(hfreq);
                        plot(hfreq,time_part,pd.freq(sel,1:ts))
                        if opts.FreqDet && ~sys.f.UseFrequencyDetector
                            plot(hfreq,time_part,pd.freq_det(sel,1:ts),'r--');
                        end
                    end

                    if opts.Ext
                        cla(hext);
                        plot(hext, time_part, pd.uneum(1:ts));
    %                     cla(hext_f);
    %                     plot(hext_f, time_part, pd.forcefactor*pd.uneum(1:ts));
                    end
                
                end
                
                if opts.Vid
                    vw.writeVideo(getframe(gcf));
                else
                    drawnow;
                end
            end
            
            if opts.Vid
                vw.close;
            end

            if isempty(opts.PM)
                pm.done;
            end
        end
    end
    
    methods(Access=protected)
        function [pd, t, y] = updatePlotData(this, pd, opts, t, y)
            mc = this.Config;
            sys = this.System;
            nf = length(mc.FibreTypes);
            
            if (opts.Freq && sys.f.UseFrequencyDetector) || opts.FreqDet
                fd = sys.f.FrequencyDetector;
                nt = length(t);
                freq = zeros(nf,nt);
                fd.reset;
                fd.WindowSize = 4;
                for i=1:nt
                    fd.processSignal(t(i),y(sys.f.moto_sarco_link_moto_out,i)');
                    freq(:,i) = fd.Frequency;
                end
                if ~isempty(opts.F)
                    freq = freq(:,1:opts.F:end);
                end
                pd.freq = freq;
                if opts.FreqDet
                    pd.freq_det = freq;
                end
            end
            
            [pd, t, y] = updatePlotData@muscle.MusclePlotter(this, pd, opts, t, y);
            nt = length(t);
            
            pd.ext_mean_current = sys.mu(4)*sys.Inputs{2,sys.inputidx}(t);
            
            if opts.Moto
                pos = sys.off_moto + (2:6:6*nf);
                pd.moto_vm = y(pos,:);
            end
            
            if opts.Sarco
                pos = sys.off_sarco + (1:56:56*nf);
                pd.sarco_pot = y(pos,:);
                
                pos = sys.off_sarco + (53:56:56*nf);
                force = bsxfun(@plus, -sys.sarco_mech_signal_offset, y(pos,:));
                force = bsxfun(@times,mc.forces_scaling,force);
                pd.sarco_force = force;
            end
            
            max_moto_signals = polyval(sys.upperlimit_poly,mc.FibreTypes);
            eff_mean_current = zeros(nf,nt);
            if opts.Spin
                pos = sys.off_spindle + (9:9:9*nf);
                pd.spindle_lambda = y(pos,:);
            
                % Freq also uses afferents if kernel expansions are used
                if opts.Aff || (opts.Freq && ~sys.f.UseFrequencyDetector)
                    afferents = zeros(2*nf,nt);
                    spindle_single_mean_current = zeros(nf,nt);
                    
                    for k=1:nf
                        spindle_pos = sys.off_spindle + (k-1)*9 + (1:9);
                        af_pos = (k-1)*2 + (1:2);
                        yspindle = y(spindle_pos,:);
                        afferents(af_pos,:) = sys.Spindle.getAfferents(yspindle);
                        spindle_single_mean_current(k,:) = sys.f.SpindleAffarentWeights*afferents(af_pos,:);
                    end
                    spindle_mean_current = mean(spindle_single_mean_current,1);
                    for k=1:nf
                        eff_mean_current(k,:) = min(max_moto_signals(k),spindle_mean_current+pd.ext_mean_current);
                    end
                    pd.afferents = afferents;
                    pd.spindle_single_mean_current = spindle_single_mean_current;
                    pd.spindle_mean_current = spindle_mean_current;
                    
                end
            else
                for k=1:nf
                    eff_mean_current(k,:) = min(max_moto_signals(k),pd.ext_mean_current);
                end
            end
            pd.eff_mean_current = eff_mean_current;
            
            if opts.Freq && ~sys.f.UseFrequencyDetector
                freq = zeros(nf,nt);
                for k=1:nf
                    x = [ones(size(pd.eff_mean_current(k,:)))*mc.FibreTypes(k);...
                        pd.eff_mean_current(k,:)];
                    freq(k,:) = sys.f.freq_kexp.evaluate(x);
%                     freq = reshape(freq,nf,[]);
                end
                pd.freq = freq;
            end
            
            % External neumann forces
            if opts.Ext
                pd.uneum = sys.Inputs{1,sys.inputidx}(t)*sys.mu(3);
%                 pd.forcefactor = norm(pd.meanforces(:,1));
            end
        end
        
        function opts = parsePlotArgs(this, args)
            %       Freq  FreqDet Aff   Geo   Moto Sarco
%             defs = [false false false false true true false];
            defs = [true true true true true true true true false];
            opts = parsePlotArgs@muscle.MusclePlotter(this, args);
            i = inputParser;
            i.KeepUnmatched = true;
            i.addParamValue('Freq',defs(1),@(v)islogical(v));
            i.addParamValue('FreqDet',defs(2),@(v)islogical(v));
            i.addParamValue('Aff',defs(3),@(v)islogical(v));
            i.addParamValue('Geo',defs(4),@(v)islogical(v));
            i.addParamValue('Moto',defs(5),@(v)islogical(v));
            i.addParamValue('Sarco',defs(6),@(v)islogical(v));
            i.addParamValue('Spin',defs(7),@(v)islogical(v));
            i.addParamValue('Ext',defs(8),@(v)islogical(v));
            i.addParamValue('GeoOnly',defs(9),@(v)islogical(v));
            i.addParamValue('MU',1:length(this.Config.FibreTypes));
            i.parse(args{:});
            opts = Utils.copyStructFields(i.Results, opts);
            if ~this.System.HasSpindle
                opts.Spin = false;
            end
        end
    end
    
end