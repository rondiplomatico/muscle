mc = GM5;

m = muscle.Model(mc);

% Experiments
experi_d = [-6:-4 -2:5]; % [mm] displacements
% Measured forces [N]
experi_F = [1.065 2.1067 3.5734 7.1639 8.8006 9.6441 10.7124 11.1749 11.295 10.8832 10.3544];
			
locations = experi_d; %[mm]
nl = length(locations);
forces = zeros(nl,length(m.Times));
forcesall = zeros(nl,length(m.Times));
pi = ProcessIndicator('Performing %d isometric tests',nl,false,nl);
for k = 1:nl
    mc.target_dist = locations(k);
    m.setConfig(mc);
    [t,y] = m.simulate;
    o = mc.getOutputOfInterest(m,t,y);
    forces(k,:) = o(1,:);
    forcesall(k,:) = o(2,:);
    pi.step;
end
pi.stop;

% Return [N] values
forces = forces/1000;
forcesall = forcesall/1000;

pm = PlotManager;%(false,1,2);
ax = pm.nextPlot('force-velocity','Forces','Displacement [mm]','Force [N]');
plot(experi_d,forces(:,end),'b');
hold(ax,'on');
plot(experi_d,experi_F,'r-x');
legend(ax,'Model','Experiment','Location','NorthWest');

% ax = pm.nextPlot('y_dir','Forces','Displacement [mm]','Force [N]');
% plot(experi_d,forcesall(:,end),'b');
% hold(ax,'on');
% plot(experi_d,experi_F,'r-x');
% legend(ax,'Model','Experiment','Location','NorthWest');
% 
% pm.done;

pm.UseFileTypeFolders = false;
pm.SaveFormats = {'eps','png'};
pm.savePlots(GM5.OutputDir);