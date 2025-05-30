% Cycle-to-cycle (C2C) JV measurement

%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - read in data files
input = 'Input_files/pmpi_v2.csv';
par = pc(input);
par.tmesh_type = 'linear'; % log10
par = refresh_device(par); % refresh the device
cycle = 1; % cycle value
xpos = 0; % define the device position point

% - - - - - - - - - - do equilibrate
soleq = equilibrate(par);

% - - - - - - - - - - do measurements
sol = doCV(soleq.ion, 0, 0, -1, 1, 1e-1, cycle, 500);

%% - - - - - - - - - - plot
dfplot_ionic.c2c(sol, xpos, cycle);
