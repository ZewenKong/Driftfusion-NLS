%% - - - - - - - - - - CODE START - - - - - - - - - -

% read in data files
input = 'Input_files/pmpi.csv';

% cv measurement setting
light_intensity = 0;
V0 = 0; V_max = -1; V_min = 1;
scan_rt = 0.2;
cycle = 1;
time_pnts = 1000;

%% (1) butler-volmer and has ionic recombination (B_ionic)
par_bvrc = pc(input);
par_bvrc.B_ionic = [0, 0, 0, 1e-14, 0];
par_bvrc.ECMswitch = 'ECMON';
% do equilibrate
soleq_bvrc = equilibrate(par_bvrc);
% transfer to simulation model
soleq_bvrc.el.par.isEquilibrate = "sim";
% soleq_bvrc.el.par.ECMswitch = "ECMON";
soleq_bvrc.ion.par.isEquilibrate = "sim";
% soleq_bvrc.ion.par. ECMswitch = "ECMON";
% do measurements
sol_bvrc = doCV(soleq_bvrc.ion, light_intensity, V0, V_max, V_min, scan_rt, cycle, time_pnts);

%% (2) no butler-volmer, but has ionic recombination
% par_rc = pc(input);
% % refresh the device
% par_rc.j0 = 0;
% par_rc = refresh_device(par_rc);
% % do equilibrate
% soleq_rc = equilibrate(par_rc);
% % transfer to simulation model
% soleq_rc.el.par.isEquilibrate = "sim";
% soleq_rc.ion.par.isEquilibrate = "sim";
% % double check
% soleq_rc.el.par.j0 = 0;
% soleq_rc.ion.par.j0 = 0;
% % do measurements
% sol_rc = doCV(soleq_rc.ion, light_intensity, V0, V_max, V_min, scan_rt, cycle, time_pnts);

%% (3) no butler-volmer, no ionic recombination
par = pc(input);
% refresh the device
par.j0 = 0;
% par.mobseti = 0;
par.mu_c = [0, 0 2.25e-14, 0, 0];
par.mu_a = [0, 0, 2.25e-14, 0, 0];
par.B_ionic = [0, 0, 0, 0, 0]; % three-layer (two-IF)
par = refresh_device(par);
% do equilibrate
soleq = equilibrate(par);
% transfer to simulation model
soleq.el.par.isEquilibrate = "sim";
soleq.ion.par.isEquilibrate = "sim";
% double check
soleq.el.par.j0 = 0;
soleq.ion.par.j0 = 0;
% do measurements
sol = doCV(soleq.ion, light_intensity, V0, V_max, V_min, scan_rt, cycle, time_pnts);

%% Plot

% plot setting
x1 = 0; x2 = 3.3e-5; % the whold device
xpos = 0;
area_coeff = 1;
use_abs = 1;
% labels = {'sol\_bv\_B_ionic', 'sol\_B_ionic', 'sol'};
% sols = cell(1, 3);
% sols{1} = sol_bvrc;
% sols{2} = sol_rc;
% sols{3} = sol;

labels = {'sol_bv', 'sol'};
sols = cell(1, 2);
sols{1} = sol_bvrc;
sols{2} = sol;

% plot
figure('Name', 'bulk perovskite device c.f.');
hold on;

for i = 1:length(sols)
    J = df_analysis.calcJ(sols{i});
    Vapp = df_analysis.calcVapp(sols{i});
    t = sols{i}.t;
    xmesh = sols{i}.x;
    ppos = getpointpos(xpos, xmesh);
    plot(Vapp, abs(J.tot(:, ppos)), 'DisplayName', labels{i});
end

hold off;
xlabel('Applied Voltage, Vapp [V]');
ylabel('Current Density, J [A cm^{-2}]');
legend('Location', 'northwest', 'FontSize', 10);
yscale log; % logrithm
