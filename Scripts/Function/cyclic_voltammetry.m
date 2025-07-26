% - - - - - - - - - - CODE START - - - - - - - - - -

%% Simulation

% parameter setting
input = 'Input_files/mapi.csv';
par = pc(input);
par.j0 = 0;
par.B_ionic = [0,0,0,0,0];
par = refresh_device(par);

% do equilibrate
soleq = equilibrate(par);

% measurements setting
light_intensity = 0;
V0 = 0;
Vmax = 0.8;
Vmin = -0.8;
scan_rt = 0.2;
cycle = 1;
time_pnts = 500;

% transfer to simulation model
soleq.el.par.isEquilibrate = "sim";
soleq.ion.par.isEquilibrate = "sim";

% do measurement
sol = doCV(soleq.ion, light_intensity, V0, Vmax, Vmin, scan_rt, cycle, time_pnts);

%% Plot

% plot setting
xpos = 0; % the examined position point
current_type = "total"; % the examined current type
area_coeff = 1; % 1 mm^2 device (0.01)
use_abs = 0; % absoulte
use_log = 0; % logarithmic

% do plot
df_plot.CVappTrack(sol, xpos, current_type, area_coeff, use_abs, use_log);
%%
df_plot.c2c(sol, xpos, area_coeff, cycle, use_abs)


