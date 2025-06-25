%% - - - - - - - - - - CODE START - - - - - - - - - -

% read in data files
input = 'Input_files/pmpi_v2.csv';
par = pc(input);
par.tmesh_type = 'log10'; % linear
par = refresh_device(par); % refresh the device

% parameters setting
cycle = 1; % cycle value
xpos = 0; % define the device position point
area_coeff = 1;
use_abs = 1;

% do equilibrate (BV)
soleq_BV = equilibrate(par);

% do measurements (BV)
sol_BV = doCV(soleq_BV.ion, 0, 0, -1, 1, 0.5, cycle, 500);

% do equilibrate
par = pc(input);
par.tmesh_type = 'log10';

par.j0 = 0; % butler-volmer
par.mobseti = 0; % ionic mobility
par = refresh_device(par);

soleq = equilibrate(par);

% do measurements
sol = doCV(soleq.ion, 0, 0, -1, 1, 0.5, cycle, 500);

% - - - - - - - - - -
sols = cell(1, 2);

sols{1} = sol_BV;
sols{2} = sol;

dfplot_ionic.d2d_var(sols, xpos);
