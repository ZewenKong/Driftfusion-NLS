%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - input file and data
input = 'Input_files/pcbm_only.csv';
par = pc(input);

xpos = 0;
cycle = 1;
area_coeff = 0.01; % compare with cm^-2
use_abs = 1; % 1 is on, 0 is off

% - - - - - - - - - - parameters setting
par.tmesh_type = 'linear';
% par.j0 = 0; % butler-volmer
% par.mobseti = 0; % ionic mobility
par = refresh_device(par);

% - - - - - - - - - - do equilibrate
soleq = equilibrate(par);

% - - - - - - - - - - do measurements
sol = doCV(soleq.ion, 0, 0, -1, 1, 0.5, 1, 500);

%% - - - - - - - - - - C2C I-V plot
dfplot_ionic.c2c(sol, xpos, area_coeff, cycle, use_abs);

%% - - - - - - - - - - Eng lvl plot
% dfplot_ionic.rELx(sol);
