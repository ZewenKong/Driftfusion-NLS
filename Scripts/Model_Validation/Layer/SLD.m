%% - - - - - - - - - - CODE START - - - - - - - - - -

% SLD, single-layer device
input = 'Input_files/pcbm_only.csv';
par = pc(input);

% parameters setting
par.tmesh_type = 'linear';
% par.j0_right = 0; 
% par.j0_left = 0;
% par.mobseti = 0;
par = refresh_device(par);

% plot setting
cycle = 1;

% do equilibrate
soleq = equilibrate(par);

%% check equilibrium solution
% dfplot.ELxnpxacx(soleq.ion);

%%
soleq.el.par.isEquilibrate = "sim";
soleq.ion.par.isEquilibrate = "sim";

%% do measurements
sol = doCV(soleq.ion, 0, 0, 1, -1, 0.5, cycle, 500); % voltage scan range adjust: 1.3, -0.7

%% do plotting

% plot setting
xpos = 4.999e-7; % position close to the right boundary
area_coeff = 1; % 1 mm^2 device (0.01)
use_abs = 0;
current_type = "total"; % current_type = "total";
absolutely = 0; % 1, turn on; 0, turn offrh
logarithmically = 0; % 1, turn on; 0, turn off

% current-voltage
dfplot_ionic.CVapp(sol, xpos, current_type, area_coeff, absolutely, logarithmically);

%%
% dfplot.JVapp(sol, xpos);
dfplot.ELxnpxacx(sol, 8);

%%
dfplot.ELxnpxacx(soleq.ion);

% 
% calculate the 'calc', by using the sol.u,x,t;
% calculate the equilibrate function