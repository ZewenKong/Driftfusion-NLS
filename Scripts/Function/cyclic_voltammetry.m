% - - - - - - - - - - CODE START - - - - - - - - - -

%% Simulation

% parameter setting
input = 'Input_files/single_layer.csv';
par = pc(input);
par_bv.j0 = 1e-3;
par_bv.RelTol = 1e-6;
par_bv.AbsTol = 1e-8;
par.B_ionic = 1e-12;
par = refresh_device(par);

% do equilibrate
soleq = equilibrate(par);

% measurements setting
light_intensity = 0;
V0 = 0;
Vmax = 1;
Vmin = -1;
scan_rt = 0.25;
cycle = 1;
time_pnts = 500;

% transfer to simulation model
soleq.el.par.isEquilibrate = "sim";
soleq.ion.par.isEquilibrate = "sim";

% do measurement
sol = doCV(soleq.ion, light_intensity, V0, Vmax, Vmin, scan_rt, cycle, time_pnts);

%% Plot

% plot setting
xpos = 3.2999e-5; % the examined position point
current_type = "total"; % the examined current type
area_coeff = 1; % 1 mm^2 device (0.01)
use_abs = 1; % absoulte
use_log = 1; % logarithmic

% do plot
df_plot.CVappTrack(sol, xpos, current_type, area_coeff, use_abs, use_log);

% df_plot.c2c(sol, xpos, area_coeff, cycle, use_abs)

% df_plot.JVapp(sol, xpos);
