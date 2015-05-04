classdef MusclePlotter < handle
% MusclePlotter: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2014-09-25
%
% @new{0,7,dw,2014-09-25} Added this class.
    
    properties
        GeoView = [46 30];
        DefaultArgs = {};
    end

    properties(SetAccess=private)
        System;
        Config;
    end
    
    properties(Access=protected)
        plotdata;
    end
    
    methods
        function this = MusclePlotter(sys)
            this.System = sys;
            this.Config = sys.Model.Config;
            this.plotdata = this.initPlotData;
        end
        
        function [pm, h] = plot(this, t, y_dofs, varargin)
            opts = this.parsePlotArgs(varargin);

            mc = this.Config;

            if isempty(opts.PM)
                if ~isempty(mc.Pool) && opts.Pool
                    pm = PlotManager(false,2,1);
                else
                    pm = PlotManager;
                end
                pm.LeaveOpen = true;
            else
                pm = opts.PM;
            end

            [pd, t] = this.updatePlotData(this.plotdata, opts, t, y_dofs);
            this.plotdata = pd;

            if opts.Vid
                [p,n,e] = fileparts(opts.Vid);
                if ~isempty(p) && exist(p,'file') ~= 7
                    error('Directory %s not found',p);
                end
                if isempty(p)
                    p = pwd;
                end
                if isempty(n)
                    n = 'output';
                end
                if isempty(e)
                    e = '.avi';
                end
                avifile = fullfile(pwd,[n e]);
                vw = VideoWriter(avifile);
                vw.FrameRate = 25;
                vw.open;
            end

            %% Loop over time
            if ~isempty(mc.Pool) && opts.Pool
                hf = pm.nextPlot('force','Activation force','t [ms]','alpha');
                axis(hf,[0 t(end) 0 1]);
                hold(hf,'on');
            end

            h = pm.nextPlot('geo',sprintf('Deformation at t=%g',t(end)),'x [mm]','y [mm]');
            zlabel(h,'z [mm]');
            axis(h, pd.geo_plotbox);
            daspect([1 1 1]);
            view(h, this.GeoView);
            hold(h,'on');

            for ts = 1:length(t)
                % Quit if figure has been closed
                if ~ishandle(h)
                    if opts.Vid
                        vw.close;
                    end
                    break;
                end
                this.plotGeometry(h, t(ts), pd.yfull(:,ts), ts, opts);

                if ~isempty(mc.Pool) && opts.Pool
                    times = t(1:ts);
                    alpha = mc.Pool.getActivation(times);
                    walpha = mc.FibreTypeWeights(1,:,1) * alpha;
                    cla(hf);
                    plot(hf,times,alpha);
                    plot(hf,times,walpha,'LineWidth',2);
                end
                
%                 % Quit if figure has been closed (extra check - runtime
%                 % thing)
%                 if ~ishandle(h)
%                     if opts.Vid
%                         vw.close;
%                     end
%                     break;
%                 end

                if ~isempty(opts.Vid) && ishandle(h)
                    vw.writeVideo(getframe(gcf));
                else
                    pause(opts.Pause);
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
        function plotGeometry(this, h, t, y_dofs, ts, opts)
            title(h,sprintf('Deformation at t=%g',t));

            pd = this.plotdata;
            sys = this.System;

            u = reshape(y_dofs(1:pd.vstart-1),3,[]);
            v = reshape(y_dofs(pd.vstart:pd.pstart-1),3,[]);
            cla(h);

            geo = this.Config.PosFE.Geometry;
            if opts.Skel
                e = pd.e;
                plot3(h,u(1,pd.no_bc),u(2,pd.no_bc),u(3,pd.no_bc),'r.','MarkerSize',14);
                for k=1:size(e,1)
                    plot3(h,u(1,[e(k,1) e(k,2)]),u(2,[e(k,1) e(k,2)]),u(3,[e(k,1) e(k,2)]),'r');
                end
            else
                p = patch('Faces',geo.PatchFaces,'Vertices',u',...
                    'Parent',h,'FaceAlpha',.3);
                if sys.HasTendons
                    set(p,'FaceVertexCData',pd.vertexcol,'FaceColor',...
                        'interp','FaceAlpha',.5,'EdgeColor','interp');
                else
                    set(p,'EdgeColor',.8*pd.musclecol,'FaceColor',pd.musclecol);
                end
            end

            % Velocities
            if opts.Velo
                quiver3(h,u(1,:),u(2,:),u(3,:),v(1,:),v(2,:),v(3,:),1,'b', 'MarkerSize',14);
            end

            %% Dirichlet conditions
            % Displacement
            u3dir = u(:,pd.bc_dir_3pos_applies);
            plot3(h,u3dir(1,:),u3dir(2,:),u3dir(3,:),'.','MarkerSize',20,'Color',[0 0 0]);
            u2dir = u(:,pd.bc_dir_2pos_applies);
            plot3(h,u2dir(1,:),u2dir(2,:),u2dir(3,:),'.','MarkerSize',20,'Color',[.5 .5 .5]);
            u1dir = u(:,pd.bc_dir_1pos_applies);
            plot3(h,u1dir(1,:),u1dir(2,:),u1dir(3,:),'.','MarkerSize',20,'Color',[.7 .7 .7]);

            % Velocity
            uvdir = u(:,pd.bool_v_bc_nodes_applies);
            plot3(h,uvdir(1,:),uvdir(2,:),uvdir(3,:),'g.','MarkerSize',20);

            %% Dirichlet Forces
            if ~isempty(pd.DF)
                residuals = pd.residuals;
                have_residuals = pd.have_residuals;
                udir = u(:,have_residuals);
                %residuals(pd.residuals_pos) = pd.DF(pd.sortidx,ts)/pd.maxdfval;
                residuals(pd.residuals_pos) = pd.DF(:,ts)/pd.maxdfval;
                quiver3(h,udir(1,:),udir(2,:),udir(3,:),...
                    residuals(1,have_residuals),residuals(2,have_residuals),...
                    residuals(3,have_residuals),1,'k', 'MarkerSize',14);
            end

            %% Neumann condition forces
            if opts.Forces
                % Plot force vectors at nodes
                uforce = u(:,pd.forces_apply);
                quiver3(h,uforce(1,:),uforce(2,:),uforce(3,:),...
                    pd.forces(1,:),pd.forces(2,:),pd.forces(3,:),0.1,'Color',[.8 .8 1]);
                % Plot force vectors at face centers
                for k=1:pd.numfaceswithforce
                    masterfacenodeidx = geo.MasterFaces(pd.force_elem_face_idx(2,k),:);
                    facenodeidx = geo.Elements(pd.force_elem_face_idx(1,k),masterfacenodeidx);

                    facecenter = mean(u(:,facenodeidx),2);
                    quiver3(h,facecenter(1),facecenter(2),facecenter(3),...
                        pd.meanforces(1,k),pd.meanforces(2,k),pd.meanforces(3,k),0.1,'LineWidth',2,'Color','b','MaxHeadSize',1);

                    if ~isempty(pd.NF)
                        pd.residual_neumann_forces(sys.bc_neum_forces_nodeidx) = pd.NF(:,ts);
                        meanforce = mean(pd.residual_neumann_forces(:,facenodeidx),2);
                        quiver3(h,facecenter(1),facecenter(2),facecenter(3),...
                        meanforce(1),meanforce(2),meanforce(3),0.1,'LineWidth',2,'Color','k','MaxHeadSize',1);
                    end
                end
            end

            %% Pressure
            if opts.Pressure
                p = y_dofs(pd.pstart:sys.num_uvp_glob);
                pn = sys.idx_p_to_u_nodes;
                pneg = p<0;
                % Negative pressures
                if any(pneg)
                    scatter3(h,u(1,pn(pneg)),u(2,pn(pneg)),u(3,pn(pneg)),...
                        -p(pneg)*10+1,[188 188 255]/255);%
                end
                % Positive pressures
                if any(~pneg)
                    scatter3(h,u(1,pn(~pneg)),u(2,pn(~pneg)),u(3,pn(~pneg)),...
                        p(~pneg)*10+1,[255 188 188]/255); %[244 149 25]/255
                end
            end

            %% a0 fibres & gp values
            fibres = sys.HasFibres && opts.Fibres;
            doI1 = any(opts.Invariants == 1);
            doLam = opts.Lambdas;
            doMR = opts.MR;
            doTMR = opts.MTRatio;
            if fibres || ~isempty(opts.Invariants) || doLam || doMR || doTMR
                for m = 1:geo.NumElements
                    u = y_dofs(1:pd.vstart-1);
                    u = u(sys.idx_u_glob_elems(:,:,m));
                    gps = u*pd.Ngp;
                    
                    if fibres
                        anull = u*sys.dNa0(:,:,m);
                        quiver3(gps(1,:),gps(2,:),gps(3,:),anull(1,:),anull(2,:),anull(3,:),.5,'.','Color','w');
                    end
                    
                    %% tendon ratio
                    if ~isempty(pd.gpcol)
                        scatter3(h,gps(1,:),gps(2,:),gps(3,:),1,pd.gpcol(:,:,m),'.','SizeData',30);
                    end
                    
                    %% Invariants
                    if doI1
                        values = sprintfc('I_1=%3.3g',pd.I1(:,m,ts));
                        text(gps(1,:),gps(2,:),gps(3,:),values,'Parent',h)
                    end
                    %% Lambdas
                    if doLam
                        %values = sprintfc('\\lambda=%3.3g',pd.lambdas(:,m,ts));
                        values = sprintfc('%3.4g',pd.lambdas(:,m,ts));
                        text(gps(1,:),gps(2,:),gps(3,:),values,'Parent',h)
                    end
                    %% Mooney-Rivlin tensors
                    if doMR
                        values = cellfun(@(m)sprintf('%g ',diag(m)), pd.MR(:,m,ts),'UniformOutput',false);
                        text(gps(1,:),gps(2,:),gps(3,:),values,'Parent',h)
                    end
                    
                    %% Tendon-Muscle ratio
                    if doTMR
                        values = sprintfc('%3.4g',sys.MuscleTendonRatioGP(:,m,ts));
                        %values = sprintfc('%3.4g',diags{:});
                        text(gps(1,:),gps(2,:),gps(3,:),values,'Parent',h)
                    end
                end
            end
        end
        
        function [pd, t, y] = updatePlotData(this, pd, opts, t, y)
            mc = this.Config;
            sys = this.System;
            dfem = mc.PosFE;
            geo = dfem.Geometry;

            pd.DF = opts.DF;
            pd.NF = opts.NF;
            %% "Speedup" factor for faster plots
            if ~isempty(opts.F)
                t = t(1:opts.F:end);
                y = y(:,1:opts.F:end);
                if ~isempty(pd.DF)
                    pd.DF = pd.DF(:,1:opts.F:end);
                end
                if ~isempty(pd.NF)
                    pd.NF = pd.NF(:,1:opts.F:end);
                end
            end
            
            %% Re-add the dirichlet nodes for geometry plotting
            pd.yfull = sys.includeDirichletValues(t, y);
            pd.geo_plotbox = this.getPlotBox(pd.yfull);

            %% Dirichlet plotting
            if ~isempty(opts.DF)
                pd.maxdfval = max(abs(opts.DF(:)))/10;
            end

            %% Forces plotting
%             if opts.Forces
                % Get forces on each node in x,y,z directions
                % This is where the connection between plane index and
                % x,y,z coordinate is "restored"
                forces = zeros(size(sys.bool_u_bc_nodes));
                forces(sys.bc_neum_forces_nodeidx) = sys.bc_neum_forces_val;

                % Forces at face centers
                force_elem_face_idx = geo.Faces(:,sys.FacesWithForce);
                pd.numfaceswithforce = size(force_elem_face_idx,2);
                meanforces = zeros(3,pd.numfaceswithforce);
                for k=1:pd.numfaceswithforce
                    masterfacenodeidx = geo.MasterFaces(force_elem_face_idx(2,k),:);
                    facenodeidx = geo.Elements(force_elem_face_idx(1,k),masterfacenodeidx);
                    meanforces(:,k) = mean(forces(:,facenodeidx),2);
                end
                pd.meanforces = meanforces;
                pd.force_elem_face_idx = force_elem_face_idx;

                % Also save forces at x,y,z nodes for plotting
                pd.forces_apply = sum(abs(forces),1) ~= 0;
                pd.forces = forces(:,pd.forces_apply);

                if ~isempty(opts.NF)
                    pd.residual_neumann_forces = zeros(size(sys.bool_u_bc_nodes));
                end
%             end

            %% Skeleton plotting
            if ~opts.Skel
        %                 light('Position',[1 1 1],'Style','infinite','Parent',h);
                pd.musclecol = [0.854688, 0.201563, 0.217188];
            end
            
            %% Show invariants or lambda stretch values
            if ~isempty(opts.Invariants) || opts.Lambdas || opts.MR
                nt = length(t);
                doI1 = any(opts.Invariants == 1);
                doLam = opts.Lambdas;
                doMR = opts.MR;
                if doI1 || doMR
                    pd.I1 = zeros(dfem.GaussPointsPerElem,geo.NumElements,nt);
                end
                if doLam
                    pd.lambdas = zeros(dfem.GaussPointsPerElem,geo.NumElements,nt);
                end
                if doMR
                    pd.MR = cell(dfem.GaussPointsPerElem,geo.NumElements,nt);
                    c10 = sys.MuscleTendonParamc10;
                    c01 = sys.MuscleTendonParamc01;
                end
                num_gp = dfem.GaussPointsPerElem;
                for m = 1:geo.NumElements
                    for gp = 1:num_gp
                        pos = 3*(gp-1)+1:3*gp;
                        dtn = dfem.transgrad(:,pos,m);        
                        elemidx_u = sys.idx_u_glob_elems(:,:,m); 
                        for ts = 1:nt
                            % Deformation gradient
                            yts = pd.yfull(:,ts) ;
                            F = yts(elemidx_u) * dtn;
                            C = F'*F;
                            if doI1 || doMR
                                pd.I1(gp,m,ts) = C(1,1)+C(2,2)+C(3,3);
                            end
                            if doLam
                                fibres = sys.a0Base(:,:,(m-1)*num_gp + gp);
                                pd.lambdas(gp,m,ts) = norm(F*fibres(:,1));
                            end
                            if doMR
                                pd.MR{gp,m,ts} = sys.MooneyRivlinICConst(gp,m)*eye(3) ...
                                    + 2*(c10(gp,m) + pd.I1(gp,m,ts)*c01(gp,m))*F ...
                                    - 2*c01(gp,m)*F*C;
                            end
                        end
                    end
                end
            end
        end
        
        function opts = parsePlotArgs(this, args)
            i = inputParser;
            i.KeepUnmatched = true;
            i.addParamValue('Vid',[],@(v)~isempty(v));
            i.addParamValue('Forces',false,@(v)islogical(v));
            i.addParamValue('Velo',false,@(v)islogical(v));
            i.addParamValue('Pressure',false,@(v)islogical(v));
            i.addParamValue('Fibres',true,@(v)islogical(v));
            i.addParamValue('Skel',false,@(v)islogical(v));
            i.addParamValue('Pool',false,@(v)islogical(v));
            i.addParamValue('PM',[],@(v)isa(v,'PlotManager'));
            i.addParamValue('DF',[]);
            i.addParamValue('NF',[]);
            i.addParamValue('F',[]);
            i.addParamValue('Invariants',[]);
            i.addParamValue('Lambdas',false);
            i.addParamValue('MR',false);
            i.addParamValue('MTRatio',false);
            i.addParamValue('Pause',.01);
            
            args = [this.DefaultArgs args];
            i.parse(args{:});
            opts = i.Results;
            if ~isempty(opts.NF)
                opts.Forces = true;
            end
        end
        
        function pd = initPlotData(this)
            sys = this.System;
            pd = struct;
            hlp = sum(sys.bool_u_bc_nodes,1);
            pd.bc_dir_3pos_applies = hlp == 3;
            pd.bc_dir_2pos_applies = hlp == 2;
            pd.bc_dir_1pos_applies = hlp == 1;
            pd.bc_dir_pos_applies = hlp >= 1; 
            pd.bool_v_bc_nodes_applies = sum(sys.bool_v_bc_nodes,1) >= 1;
            pd.no_bc = ~pd.bc_dir_pos_applies & ~pd.bool_v_bc_nodes_applies;

            mc = this.Config;
            dfem = mc.PosFE;
            geo = dfem.Geometry;
            posdofs = geo.NumNodes * 3;
            pd.vstart = posdofs+1;
            pd.pstart = 2*posdofs+1;
            pd.e = geo.Edges;

            %% Dirichlet forces
            pd.have_residuals = pd.bc_dir_pos_applies | pd.bool_v_bc_nodes_applies;
            % By sorting and combining the pos/velo Dir BC, the
            % plotting of mixed BCs on one node is plotted correctly.
            pd.residuals_pos = sys.bool_u_bc_nodes | sys.bool_v_bc_nodes;
            [~, pd.sortidx] = sort([sys.idx_u_bc_glob; sys.idx_v_bc_glob-posdofs]);
            % Preallocate the residuals matrix
            pd.residuals = zeros(size(pd.residuals_pos));

            %% Fibres
            %if sys.HasFibres || ~isempty(opts.Invariants)
                pd.Ngp = dfem.N(dfem.GaussPoints);
            %end
            
            % Set muscle color either way
            pd.musclecol = [0.854688, 0.201563, 0.217188];
            pd.gpcol = [];
            %% Muscle/Tendon patch vertex colors
            if sys.HasTendons
                tendoncol = [1, .8, .6]; %Micha [.9 .7 .5]
                % Vertex coloring
                pd.vertexcol = ones(geo.NumNodes,1)*pd.musclecol ...
                    + sys.MuscleTendonRatioNodes'*(tendoncol-pd.musclecol);
                % Gauss point coloring
                tmr = sys.MuscleTendonRatioGP;
                pd.gpcol(:,:,1) = pd.musclecol(1) + tmr*(tendoncol(1)-pd.musclecol(1));
                pd.gpcol(:,:,2) = pd.musclecol(2) + tmr*(tendoncol(2)-pd.musclecol(2));
                pd.gpcol(:,:,3) = pd.musclecol(3) + tmr*(tendoncol(3)-pd.musclecol(3));
                pd.gpcol = permute(pd.gpcol,[1 3 2]);%*.8;
            end
        end
%     end
%     
%     methods(Access=private)
        function box = getPlotBox(this, uvw)
            xpos = 1:3:this.Config.PosFE.Geometry.NumNodes*3;
            box = [min(min(uvw(xpos,:))) max(max(uvw(xpos,:)))...
                min(min(uvw(xpos+1,:))) max(max(uvw(xpos+1,:)))...
                min(min(uvw(xpos+2,:))) max(max(uvw(xpos+2,:)))];
            diam = diff(box)*.1;
            box = box + [-diam(1) diam(1) -diam(3) diam(3) -diam(5) diam(5)];
        end
    end
    
end