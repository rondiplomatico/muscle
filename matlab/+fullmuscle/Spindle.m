classdef Spindle < KerMorObject
% Spindle: 
%
% @docupdate
%
% @author Daniel Wirtz @date 2014-09-18
%
% @new{0,7,dw,2014-09-18} Added this class.
%
% This class is part of the framework
% KerMor - Model Order Reduction using Kernels:
% - \c Homepage http://www.agh.ians.uni-stuttgart.de/research/software/kermor.html
% - \c Documentation http://www.agh.ians.uni-stuttgart.de/documentation/kermor/
% - \c License @ref licensing
    
    properties(SetAccess=private)
        spindleConst;
        JSparsityPattern;
        y0;
    end
    
    methods
        
        function this = Spindle
            this.initSpindleConst;
           
            % Sparsity pattern
            i = [1 3 4 5 2 4 3 6 4 7 5 8 3 4 5 3 4 5];
            j = [1 1 1 1 2 2 3 3 4 4 5 5 6 7 8 9 9 9 ];
            this.JSparsityPattern = sparse(i,j,true,9,9);
            
            c = this.spindleConst;
            
            % Initial values
            y = zeros(9,1);
            y(1) = 0; %0.5*0.735;
            y(2) = 0; %0.5*0.735;            
            y(3:5) = 0.000634066218078;
%             y(3) = 0.000634066218078;
%             y(4) = 0.000634067218078;
%             y(5) = 0.000634068218078;
%             % DONT use exactly identical initial values. This makes sharp
%             % comparisons between afferents insecure. Instead, slightly
%             % changed initial conditions yield clear, intuitive decisions.
            y(6) = c(1)*(c(12)-c(17)); %0.020731324015575;
            y(7) = c(31)*(c(42)-c(47));%0.020731334015575;
            y(8) = c(61)*(c(72)-c(77));%0.020731344015575;
%             y(6:8) = 0.020731324015575;
            y(9) = 0.95; % length
%             y(10) = 0; % output phase value, 'primary_afferent'
%             y(11) = 0; % output phase value, 'secondary_afferent'
            this.y0 = y;
        end
        
        % c(20) => c(50) in f_stat_bag2
        
        function dy = dydt(this, y, t, moto_sig, Ldot, Lddot)
            % dydt_spindle: Function for change rate of spindle y.

            %% Init
            moto_freq_dyn = moto_sig(1);
            moto_freq_static = moto_sig(2);
            
            C = 0.42;
            C = C + (Ldot > 0)*.58; % C = 1 if Ldot > 0
            
            c = this.spindleConst;
            dy = zeros(size(y));
            
            %% Eval
            
            tmp1 = Ldot-y(3,:)./c(1);
            tmp1neg = tmp1 < 0;
            tmp1pos = ~tmp1neg;
            tmp2 = Ldot-y(4,:)./c(31);
            tmp2neg = tmp2 < 0;
            tmp2pos = ~tmp2neg;
            tmpc = Ldot-y(5,:)./c(61);
            tmpcneg = tmpc < 0;
            tmpcpos = ~tmpcneg;
            
            % f_dyn_bag1
            dy(1,:) = (moto_freq_dyn^c(22)./(moto_freq_dyn^c(22)+c(21)^c(22)) - y(1,:))/c(20);
            % f_stat_bag2
            dy(2,:) = (moto_freq_static.^c(52)/(moto_freq_static^c(52)+c(51)^c(52)) - y(2,:))/c(50);
            % f_stat_chain
            f_stat_chain = moto_freq_static^c(82)/(moto_freq_static^c(82)+c(81)^c(82));

            %% Bag1
            b1  = c(4)+c(5)*y(1,:);
            g1 = c(7)*y(1,:);
            % y3 = T' => dy3 = T''
            dy(3,:) = c(1)/c(3)*(...
                C.*b1.*(tmp1pos.*tmp1.^c(15) - tmp1neg.*((-tmp1).^c(15)))...
                .*(y(9,:)-c(17)-y(6,:)/c(1)-c(16))...
                + c(2)*(y(9,:)-c(17)-y(6,:)/c(1)-c(18))...
                + c(3)*Lddot + g1 - y(6,:));
            % y6 = T => dy6 = T' = y3
            dy(6,:) = y(3,:);

            %% Bag2
            b2  = c(34)+c(35)*y(1,:)+c(36)*y(2,:);
            g2 = c(37)*y(1,:)+c(38)*y(2,:);
            %C.*b2.*sgntmp2.*abs(Ldot-y(4,:)/c(31)).^c(45)...
            % y4 = T' => dy4 = T''
            dy(4,:) = c(31)/c(33)*(...
                C.*b2.*(tmp2pos.*tmp2.^c(45) - tmp2neg.*((-tmp2).^c(45)))...
                .*(y(9,:) - c(47)-y(7,:)/c(31)-c(46))...
                + c(32)*(y(9,:)-c(47)-y(7,:)/c(31)-c(48))...
                + c(33)*Lddot + g2 - y(7,:));
            % y7 = T => dy7 = T' = y4
            dy(7,:) = y(4,:);

            %% Chain
            bc  = c(64)+c(65)*y(1,:) + c(66)*f_stat_chain;
            gc = c(67)*y(1,:) + c(68)*f_stat_chain;
            % y4 = T' => dy4 = T''
            dy(5,:)  = c(61)/c(63)*(...
                C.*bc.*(tmpcpos.*tmpc.^c(75) - tmpcneg.*((-tmpc).^c(75))) .* ...
                (y(9,:)-c(77)-y(8,:)/c(61)-c(76)) +c(62) .* ...
                (y(9,:)-c(77)-y(8,:)/c(61)-c(78)) +c(63) .* Lddot + gc - y(8,:));
            % y8 = T => dy8 = T' = y5
            dy(8,:)  = y(5,:);

            % dy9 = L' (y9 = L) in [L_0]
            dy(9,:) = Ldot;
        end
        
        function af  = getAfferents(this, y)
            af = zeros(2,size(y,2));
            c = this.spindleConst;
            
            % algebraic frequency from bag1,bag2 and chain
            Aff_pot_prim_bag1  = c(14)*(y(6,:)/c(1) - (c(12)-c(17)));
            Aff_pot_prim_bag2  = c(44)*(y(7,:)/c(31) - (c(42)-c(47)));
            Aff_pot_prim_chain = c(74)*(y(8,:)/c(61) - (c(72)-c(77)));
            bag2_and_chain = Aff_pot_prim_bag2 + Aff_pot_prim_chain;
            
            % Primary_afferent
            af(1,:) = max(Aff_pot_prim_bag1,bag2_and_chain) ...
                + 0.156*min(Aff_pot_prim_bag1,bag2_and_chain);
            
            Aff_pot_sec_bag2  = c(53)*(c(41)*c(49)/c(47)*(y(7,:)/c(31)-(c(42)-c(47)))...
                + (1-c(41))*c(49)/c(48)*(y(9,:)-y(7,:)/c(31)-c(47)-c(43)));
            Aff_pot_sec_chain = c(83)*(c(71)*c(79)/c(77)*(y(8,:)/c(61)-(c(72)-c(77)))...
                + (1-c(71))*c(79)/c(78)*(y(9,:)-y(8,:)/c(61)-c(77)-c(73)));
            
            % Secondary_afferent
            af(2,:) = Aff_pot_sec_bag2 + Aff_pot_sec_chain;
        end
        
        function plot(this, t, y)
            y = y';
            pm = PlotManager(false,2,2);
            pm.LeaveOpen = true;
            
            ax = pm.nextPlot('freq','Frequencies','t','f');
            plot(ax,t,y(1:2,:));
            
            ax = pm.nextPlot('length','Fibre length','t','L [L_0]');
            plot(ax,t,y(9,:));
            
            af = this.getAfferents(y);
            
            ax = pm.nextPlot('prim','Primary afferent','t','primary');
            plot(ax,t,af(1,:));
            
            ax = pm.nextPlot('second','Secondary afferent','t','secondary');
            plot(ax,t,af(2,:));
            
            pm.done;
        end
        
        function J = Jdydt(this, y, t, moto_sig, Ldot, Lddot)
            % computes the Jacobian for dydt
            %
            moto_freq_dyn = moto_sig(1);
            moto_freq_static = moto_sig(2);

            C = 0.42;
            C = C + (Ldot > 0)*.58; % C = 1 if Ldot > 0
            
            c = this.spindleConst;
            
            tmp1 = Ldot-y(3,:)./c(1);
            tmp1neg = tmp1 < 0;
            tmp1pos = ~tmp1neg;
            tmp2 = Ldot-y(4,:)./c(31);
            tmp2neg = tmp2 < 0;
            tmp2pos = ~tmp2neg;
            tmpc = Ldot-y(5,:)./c(61);
            tmpcneg = tmpc < 0;
            tmpcpos = ~tmpcneg;
                        
            J = sparse(9,9);
            J(1,1) = -1/c(20);
            J(2,2) = -1/c(50);
            J(3,1) = (c(1)*(c(7) + C*c(5)*(tmp1neg*(y(3)/c(1) - Ldot)^c(15) - tmp1pos*(Ldot - y(3)/c(1))^c(15))*(c(16) + c(17) - y(9) + y(6)/c(1))))/c(3);
            J(3,3) = (C*c(1)*((c(15)*tmp1pos*(Ldot - y(3)/c(1))^(c(15) - 1))/c(1) + (c(15)*tmp1neg*(y(3)/c(1) - Ldot)^(c(15) - 1))/c(1))*(c(4) + c(5)*y(1))*(c(16) + c(17) - y(9) + y(6)/c(1)))/c(3);
            J(3,6) = -(c(1)*(c(2)/c(1) - (C*(tmp1neg*(y(3)/c(1) - Ldot)^c(15) - tmp1pos*(Ldot - y(3)/c(1))^c(15))*(c(4) + c(5)*y(1)))/c(1) + 1))/c(3);
            J(3,9) = (c(1)*(c(2) - C*(tmp1neg*(y(3)/c(1) - Ldot)^c(15) - tmp1pos*(Ldot - y(3)/c(1))^c(15))*(c(4) + c(5)*y(1))))/c(3);
            J(4,1) = (c(31)*(c(37) + C*c(35)*(tmp2neg*(y(4)/c(31) - Ldot)^c(45) - tmp2pos*(Ldot - y(4)/c(31))^c(45))*(c(46) + c(47) - y(9) + y(7)/c(31))))/c(33);
            J(4,2) = (c(31)*(c(38) + C*c(36)*(tmp2neg*(y(4)/c(31) - Ldot)^c(45) - tmp2pos*(Ldot - y(4)/c(31))^c(45))*(c(46) + c(47) - y(9) + y(7)/c(31))))/c(33);
            J(4,4) = (C*c(31)*((c(45)*tmp2pos*(Ldot - y(4)/c(31))^(c(45) - 1))/c(31) + (c(45)*tmp2neg*(y(4)/c(31) - Ldot)^(c(45) - 1))/c(31))*(c(34) + c(35)*y(1) + c(36)*y(2))*(c(46) + c(47) - y(9) + y(7)/c(31)))/c(33);
            J(4,7) = -(c(31)*(c(32)/c(31) - (C*(tmp2neg*(y(4)/c(31) - Ldot)^c(45) - tmp2pos*(Ldot - y(4)/c(31))^c(45))*(c(34) + c(35)*y(1) + c(36)*y(2)))/c(31) + 1))/c(33);
            J(4,9) = (c(31)*(c(32) - C*(tmp2neg*(y(4)/c(31) - Ldot)^c(45) - tmp2pos*(Ldot - y(4)/c(31))^c(45))*(c(34) + c(35)*y(1) + c(36)*y(2))))/c(33);
            J(5,1) = (c(61)*(c(67) + C*c(65)*(tmpcneg*(y(5)/c(61) - Ldot)^c(75) - tmpcpos*(Ldot - y(5)/c(61))^c(75))*(c(76) + c(77) - y(9) + y(8)/c(61))))/c(63);
            J(5,5) = (C*c(61)*((c(75)*tmpcpos*(Ldot - y(5)/c(61))^(c(75) - 1))/c(61) + (c(75)*tmpcneg*(y(5)/c(61) - Ldot)^(c(75) - 1))/c(61))*(c(64) + c(65)*y(1) + (c(66)*moto_freq_static^c(82))/(c(81)^c(82) + moto_freq_static^c(82)))*(c(76) + c(77) - y(9) + y(8)/c(61)))/c(63);
            J(5,8) = -(c(61)*(c(62)/c(61) - (C*(tmpcneg*(y(5)/c(61) - Ldot)^c(75) - tmpcpos*(Ldot - y(5)/c(61))^c(75))*(c(64) + c(65)*y(1) + (c(66)*moto_freq_static^c(82))/(c(81)^c(82) + moto_freq_static^c(82))))/c(61) + 1))/c(63);
            J(5,9) = (c(61)*(c(62) - C*(tmpcneg*(y(5)/c(61) - Ldot)^c(75) - tmpcpos*(Ldot - y(5)/c(61))^c(75))*(c(64) + c(65)*y(1) + (c(66)*moto_freq_static^c(82))/(c(81)^c(82) + moto_freq_static^c(82)))))/c(63);
            J(6,3) = 1;
            J(7,4) = 1;
            J(8,5) = 1;
        end
    end
    
    methods(Access=private)
        function initSpindleConst(this)
            % initSpindleConst: Init function for spindle constants
            
            % copied from spindle_whole_2012_10_11.m, line 217 - 291
            c = zeros(1,83);

            % Non-existent for Bag_1: c(6), c(8), c(11), c(13), c(19)
            % c(9) replaced with C as an input parameter
            c(1:11) = [10.4649, 0.15, 0.0002, 0.0605, 0.2592, 0, 0.0289, 0, 0, 0, 0];
            c(12:22) = [0.0423, 0, 20000, 0.3, 0.46, 0.04, 0.76, 0, 0.149, 60, 2];

            % Non-existent for Bag_2: c(35), c(37); c(39) replaced with C as an input parameter
            c(31:42) = [10.4649, 0.15, 0.0002, 0.0822, 0, -0.046, 0, 0.0636, 0, 0, 0.7, 0.0423];
            c(43:53) = [0.89, 10000, 0.3, 0.46, 0.04, 0.76, 0.04, 0.205, 60, 2, 7250];

            % Non-existent for Chain: c(65), c(67),c(80); c(69) replaced with C as an input parameter
            c(61:72) = [10.4649, 0.15, 0.0002, 0.0822, 0, -0.069, 0, 0.0954, 0, 0, 0.7, 0.0423];
            c(73:83) = [0.89, 10000, 0.3, 0.46, 0.04, 0.76, 0.04, 0, 90, 2, 7250];

            this.spindleConst = c;
        end
    end
    
    methods(Static)
        function test_Spindle
            s = fullmuscle.Spindle;
            t = 0:.01:3.3;
            Lddot = 0;
            startt = 1;
            
%             ldot = [.11 .66 1.55]; [L_0/s]
            ldot = .11;
            lendiff = 1.08 - .95;
            moto_sig = [30 0];
            opt = odeset;
            for k = 1:length(ldot)
                endt = startt + lendiff/ldot(k);

                % opt = odeset('Jacobian',@odejac);
                % [t,y] = ode15s(@(t,y)s.dydt(y,t,moto_sig,Ldot(t),Lddot),t,s.y0,opt);

                dy0 = s.dydt(s.y0,0,moto_sig,Ldot(1),Lddot);
                opt = odeset(opt,'Jacobian',@ode15ijac);
                [t,y] = ode15i(@(t,y,dy)dy-s.dydt(y,t,moto_sig,Ldot(t),Lddot),t,s.y0,dy0,opt);

            end
            
            function [J,M] = ode15ijac(t,y,~)
                J = -s.Jdydt(y,t,moto_sig,Ldot(t),Lddot);
                M = speye(9);
            end
            
            function J = odejac(t,y)
                J = s.Jdydt(y,t,moto_sig,Ldot(t),Lddot);
            end
            
            function ld = Ldot(t)
                ld = 0;
                if t >= startt && t < endt
                    ld = ldot(k);
                end
            end

            s.plot(t,y);
        end
        
        function res = test_Spindle_Jac
            s = fullmuscle.Spindle;
            d = 9;
            n = 10;
            moto_sig = rand(n,2);
            Ldot = 0;
            Lddot = 0;
            testx = [s.y0 bsxfun(@times,s.y0,rand(d,n-1)*10)];
            res = true;
            for k = 1:n
                x = testx(:,k);
                X = repmat(x,1,d);
                dx = ones(d,1)*sqrt(eps(class(x))).*max(abs(x),1);
                DX = sparse(1:d,1:d,dx,d,d);

                % Evaluate makes use of built-in multi-argument evaluation
                Jc = (s.dydt(X+DX,0,moto_sig(k,:),Ldot,Lddot)...
                    - repmat(s.dydt(x,0,moto_sig(k,:),Ldot,Lddot),1,d))*diag(1./dx);
                J = s.Jdydt(x,0,moto_sig(k,:),Ldot,Lddot);
                absdiff = abs(Jc-J);
                reldiff = abs(absdiff./J);
                maxabsdiff = max(absdiff(s.JSparsityPattern));
                maxreldiff = max(reldiff(s.JSparsityPattern));
                fprintf('Max abs diff: %g, max rel diff: %g\n',maxabsdiff,maxreldiff);
                res = res & maxreldiff < .005;
            end
        end
        
        function createSpindleJac
            c = sym('c',[83 1]);
            y = sym('y',[11 1]); 
            dy = sym('dy',[11 1]);

            moto_freq_dyn = sym('moto_freq_dyn','positive');
            moto_freq_static = sym('moto_freq_static','positive');
            C = sym('C','positive');
            Ldot = sym('Ldot');
            Lddot = sym('Lddot');
            
            tmp1neg = sym('tmp1neg','positive');
            tmp1pos = sym('tmp1pos','positive');
            tmp2neg = sym('tmp2neg','positive');
            tmp2pos = sym('tmp2pos','positive');
            tmpcneg = sym('tmpcneg','positive');
            tmpcpos = sym('tmpcpos','positive');

            %% Definitions (copy code here as of ""Eval"" section in dydt)
            tmp1 = Ldot-y(3,:)./c(1);
%             tmp1neg = tmp1 < 0;
%             tmp1pos = ~tmp1neg;
            tmp2 = Ldot-y(4,:)./c(31);
%             tmp2neg = tmp2 < 0;
%             tmp2pos = ~tmp2neg;
            tmpc = Ldot-y(5,:)./c(61);
%             tmpcneg = tmpc < 0;
%             tmpcpos = ~tmpcneg;
            
            % f_dyn_bag1
            dy(1,:) = (moto_freq_dyn^c(22)./(moto_freq_dyn^c(22)+c(21)^c(22)) - y(1,:))/c(20);
            % f_stat_bag2
            dy(2,:) = (moto_freq_static.^c(52)/(moto_freq_static^c(52)+c(51)^c(52)) - y(2,:))/c(50);
            % f_stat_chain
            f_stat_chain = moto_freq_static^c(82)/(moto_freq_static^c(82)+c(81)^c(82));

            %% Bag1
            b1  = c(4)+c(5)*y(1,:);
            g1 = c(7)*y(1,:);
            % y3 = T' => dy3 = T''
            dy(3,:) = c(1)/c(3)*(...
                C.*b1.*(tmp1pos.*tmp1.^c(15) - tmp1neg.*((-tmp1).^c(15)))...
                .*(y(9,:)-c(17)-y(6,:)/c(1)-c(16))...
                + c(2)*(y(9,:)-c(17)-y(6,:)/c(1)-c(18))...
                + c(3)*Lddot + g1 - y(6,:));
            % y6 = T => dy6 = T' = y3
            dy(6,:) = y(3,:);

            %% Bag2
            b2  = c(34)+c(35)*y(1,:)+c(36)*y(2,:);
            g2 = c(37)*y(1,:)+c(38)*y(2,:);
            %C.*b2.*sgntmp2.*abs(Ldot-y(4,:)/c(31)).^c(45)...
            % y4 = T' => dy4 = T''
            dy(4,:) = c(31)/c(33)*(...
                C.*b2.*(tmp2pos.*tmp2.^c(45) - tmp2neg.*((-tmp2).^c(45)))...
                .*(y(9,:) - c(47)-y(7,:)/c(31)-c(46))...
                + c(32)*(y(9,:)-c(47)-y(7,:)/c(31)-c(48))...
                + c(33)*Lddot + g2 - y(7,:));
            % y7 = T => dy7 = T' = y4
            dy(7,:) = y(4,:);

            %% Chain
            bc  = c(64)+c(65)*y(1,:) + c(66)*f_stat_chain;
            gc = c(67)*y(1,:) + c(68)*f_stat_chain;
            % y4 = T' => dy4 = T''
            dy(5,:)  = c(61)/c(63)*(...
                C.*bc.*(tmpcpos.*tmpc.^c(75) - tmpcneg.*((-tmpc).^c(75))) .* ...
                (y(9,:)-c(77)-y(8,:)/c(61)-c(76)) +c(62) .* ...
                (y(9,:)-c(77)-y(8,:)/c(61)-c(78)) +c(63) .* Lddot + gc - y(8,:));
            % y8 = T => dy8 = T' = y5
            dy(8,:)  = y(5,:);

            % dy9 = L' (y9 = L) in [L_0]
            dy(9,:) = Ldot;

            %% Create partial derivatives
            JP = sparse(false(9,9));
            f = 1;
            fprintf(f,'J = sparse(9,9);\n');
            for i = 1:11
                for j = 1:11
                    pd = diff(dy(i),y(j));
                    if pd ~= 0
                        JP(i,j) = true;
                        % Convert y14 to y(14) etc, also for c
                        body = regexprep(char(pd), '(c|y)(\d+)', '$1($2)');
                        fprintf(f,'J(%d,%d) = %s;\n',i,j,body);
                    end
                end
            end
            [i,j] = find(JP);
            sprintf('%d ',i)
            sprintf('%d ',j)
        end
    end
    
end