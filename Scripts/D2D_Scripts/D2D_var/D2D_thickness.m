% Device-to-device (D2D) JV measurements,

%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files/pmpi.csv';
active_thickness = [2.2e-5, 2.7e-5, 3.2e-5];
var = active_thickness;

% - - - - - - - - - - handles
soleqs = cell(1, length(var));
sols = cell(size(soleqs));
on_off_ratio_handle = cell(size(sols));

% - - - - - - - - - - data processing
par = pc(input);
par.prob_distro_function = 'Boltz';
par.tmesh_type = 'linear';
par = refresh_device(par);
xpos = 0;

for i = 1:length(var)
    par.d(3) = var(i); % active layer thickness
    par = refresh_device(par);
    soleqs{i} = equilibrate(par);
end

% - - - - - - - - - - do measurements
for i = 1:length(soleqs)
    soleq = soleqs{i};
    sol = doCV(soleq.ion, 0, 0, -1, 1, 1e-1, 1, 500); % solution
    sols{i} = sol;
end

%% - - - - - - - - - - plot
dfplot_ionic.d2d_var(sols, xpos);
