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
            nf = length(mc.FibreTypes);
            
            if isempty(opts.PM)
                pm = PlotManager(false,2,3);
                pm.LeaveOpen = true;
            else
                pm = opts.PM;
            end
            
            [t, y] = this.updatePlotData(opts, t, y);
            
            if opts.Vid
                avifile = fullfile(pwd,'output.avi');
                vw = VideoWriter(avifile);
                vw.FrameRate = 25;
                vw.open;
            end
            
            %% Loop over time
            if ~isempty(mc.Pool) && r.Pool
                hf = pm.nextPlot('force','Activation force','t [ms]','alpha');
                axis(hf,[0 t(end) 0 1]);
                hold(hf,'on');
            end
            
%             pmgeo = PlotManager;
%             pmgeo.LeaveOpen = true;
            pmgeo = pm;
            h_geo = pmgeo.nextPlot('geo',sprintf('Deformation at t=%g',t(end)),'x [mm]','y [mm]');
            zlabel(h_geo,'z [mm]');
            axis(h_geo, this.getPlotBox(y));
            daspect([1 1 1]);
            view(h_geo, [46 30]);
            hold(h_geo,'on');
            
            h1 = pm.nextPlot('signal','Motoneuron signal','t [ms]','V_m');
            pos = sys.off_moto_full + (2:6:6*nf);
            vals = y(pos,:);
            axis(h1,[0 t(end) min(vals(:)) max(vals(:))]);
            hold(h1,'on');

            h2 = pm.nextPlot('force','Action potential','t [ms]','V_m');
            pos = sys.off_sarco_full + (1:56:56*nf);
            vals = y(pos,:);
            axis(h2,[0 t(end) min(vals(:)) max(vals(:))]);
            hold(h2,'on');

            h3 = pm.nextPlot('force','Activation','t [ms]','A_s');
            pos = sys.off_sarco_full + (53:56:56*nf);
            vals = y(pos,:);
            vals = bsxfun(@times, min(1,t), vals);
            axis(h3,[0 t(end) min(vals(:)) max(vals(:))]);
            hold(h3,'on');

            h4 = pm.nextPlot('spindle','Afferents','t [ms]','aff');
            pos = sys.off_spindle_full + (1:9*nf);
            vals = sys.Spindle.getAfferents(y(pos,:));
            axis(h4,[0 t(end) min(vals(:)) max(vals(:))]);
            hold(h4,'on');

            h5 = pm.nextPlot('spindle','Signal','t [ms]','aff');
            maxcurrents = polyval(sys.upperlimit_poly,mc.FibreTypes);
            sp_sig = sys.f.SpindleAffarentWeights*vals;
            sp_sig = min(maxcurrents(1) - sys.mu(4),sp_sig);
            axis(h5,[0 t(end) min(sp_sig(:)) max(sp_sig(:))]);
            hold(h5,'on');
            
            for ts = 1:length(t)
                % Quit if figure has been closed
                if ~ishandle(h_geo)
                    break;
                end
                
                %% Plot geometry into h_geo
                this.plotGeometry(h_geo, t(ts), y(:,ts), ts, opts);
                
                time_part = t(1:ts);
                pos = sys.off_moto_full + (2:6:6*nf);
                signal = y(pos,1:ts);
                %walpha = mc.FibreTypeWeights(1,:,1) * signal;
%                 cla(h(1));
                plot(h1,time_part,signal);
                %plot(h,times,walpha,'LineWidth',2);
                
                pos = sys.off_sarco_full + (1:56:56*nf);
                force = y(pos,1:ts);
                %walpha = mc.FibreTypeWeights(1,:,1) * signal;
%                 cla(h(2));
                plot(h2,time_part,force);
                
                pos = sys.off_sarco_full + (53:56:56*nf);
                force = y(pos,1:ts);
                force = bsxfun(@plus, -sys.sarco_mech_signal_offset, force);
                %force = bsxfun(@times, min(1,time_part), force);
                
                walpha = mc.FibreTypeWeights(1,:,1) * force;
%                force = [force; walpha]
%                 cla(h(3));
                plot(h3,time_part,force,'r',time_part,walpha,'b');
%                 plotyy(h(3),time_part,force,time_part,walpha);

                pos = sys.off_spindle_full + (1:9*nf);
                affarents = sys.Spindle.getAfferents(y(pos,1:ts));
                plot(h4,time_part,affarents');
                
                maxcurrents = polyval(sys.upperlimit_poly,mc.FibreTypes);
                sp_sig = sys.f.SpindleAffarentWeights * affarents;
                sp_sig = min(maxcurrents(1) - sys.mu(4),sp_sig);
                plot(h5,time_part,sp_sig);
                
                drawnow;
            end
            
            if opts.Vid
                vw.close;
            end

            if isempty(opts.PM)
                pm.done;
            end
        end
    end
    
end