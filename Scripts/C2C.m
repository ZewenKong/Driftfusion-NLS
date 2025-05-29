% Cycle-to-cycle (C2C) JV measurement,
% to define the endurance and the stochasticity

% - - - - - - - - - - initialisation
equilibrate_init();
initialise_df;

%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - read in data files
input = 'Input_files/pmpi_v2.csv';

par = pc(input); % original parameters
par.tmesh_type = 'linear'; % log10

par = refresh_device(par); % refresh the device

soleq = equilibrate(par);
cycle = 1; % cycle value

% - - - - - - - - - - do measurements
sol = doCV(soleq.ion, 0, 0, -1, 1, 1e-1, cycle, 500);

%% - - - - - - - - - -
J = dfana.calcJ(sol);
Vapp = dfana.calcVapp(sol);
t = sol.t;

xmesh = sol.x;
xpos = 0; % define the device position point
ppos = getpointpos(xpos, xmesh);

dfplot_ionic.c2c(J, Vapp, ppos, cycle); % plot
