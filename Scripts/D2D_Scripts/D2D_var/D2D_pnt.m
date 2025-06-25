%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files/pmpi_v2.csv';
pnts = [100, 250, 300, 400, 500, 600, 750, 800];

% - - - - - - - - - - handle
sols = cell(1, length(scan_rts)); % solutions cell array

% - - - - - - - - - - data processing
par = pc(input);
par.prob_distro_function = 'Boltz';
par.tmesh_type = 'linear';
xpos = 0;
soleq = equilibrate(par);

% - - - - - - - - - - do measurements
for i = 1:length(scan_rts)
    sol = doCV(soleq.ion, 0, 0, -1, 1, 0.5, 1, pnts(i)); % solution
    sols{i} = sol;
end

%% - - - - - - - - - - plot
dfplot_ionic.d2d_var(sols, xpos);
