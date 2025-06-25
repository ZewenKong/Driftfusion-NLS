% Cycle-to-cycle (C2C) JV measurement

%% - - - - - - - - - - CODE START - - - - - - - - - -

% read in data files
input = 'Input_files/pmpi_v2.csv';
par = pc(input);
par.tmesh_type = 'log10'; % linear
par = refresh_device(par); % refresh the device

cycle = 1; % cycle value
xpos = 0; % define the device position point
area_coeff = 1;
use_abs = 1;

% - - - - - - - - - - do equilibrate
soleq = equilibrate(par);

% - - - - - - - - - - do measurements
sol = doCV(soleq.ion, 0, 0, -1, 1, 0.5, cycle, 500);

%% - - - - - - - - - - plot
dfplot_ionic.c2c(sol, xpos, area_coeff, cycle, use_abs);
dfplot.ELx(sol);
dfplot_ionic.npxacx(sol);

