%% - - - - - - - - - - CODE START - - - - - - - - - -

% read in data files
input = 'Input_files/mapi.csv';
% cv measurement setting
light_intensity = 0;
V0 = 0; V_max = 1; V_min = -1;
scan_rt = 0.2;
cycle = 1;
time_pnts = 500;

%% ECM
par_ecm = pc(input);

par_ecm.isECM = "ecm_on";
par_ecm.k0_trap = 1e-3;
par_ecm.dynamic_adp = 1e-9;

par_ecm.j0 = 1e-3;
par_ecm.B_ionic = [0, 0, 0, 1e-12, 0];

par_ecm.RelTol = 1e-6;
par_ecm.AbsTol = 1e-8;
par_ecm = refresh_device(par_ecm);

soleq_ecm = equilibrate(par_ecm);
soleq_ecm.el.par.isEquilibrate = "sim";
soleq_ecm.ion.par.isEquilibrate = "sim";

sol_ecm = doCV(soleq_ecm.ion, light_intensity, V0, V_max, V_min, scan_rt, cycle, time_pnts);

%% General
par = pc(input);

par.mu_c = [0, 0 2.25e-14, 0, 0];
par.mu_a = [0, 0, 2.25e-14, 0, 0];
par.k0_trap = 1e-3;
par.dynamic_adp = 1e-9;

par.j0 = 0;
par.B_ionic = [0, 0, 0, 0, 0];

par.RelTol = 1e-6;
par.AbsTol = 1e-8;
par = refresh_device(par);

soleq = equilibrate(par);
soleq.el.par.isEquilibrate = "sim";
soleq.ion.par.isEquilibrate = "sim";

sol = doCV(soleq.ion, light_intensity, V0, V_max, V_min, scan_rt, cycle, time_pnts);

%% Plot
x1 = 0; x2 = 3.3e-5; % the whold device
xpos = 0;
area_coeff = 1;
use_abs = 1;

labels = {'sol-ecm', 'sol'};
sols = cell(1, 2);
sols{1} = sol_ecm;
sols{2} = sol;

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
