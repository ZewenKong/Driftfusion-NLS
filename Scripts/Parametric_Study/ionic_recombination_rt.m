%% - - - - - - - - - - CODE START - - - - - - - - - -

% data inputs
input = 'Input_files/pmpi_v2.csv';
B_ionic_array = [1e-10, 1e-12, 1e-14, 1e-16];
var = B_ionic_array;

% handles
soleqs = cell(1, length(var));
sols = cell(1, length(soleqs));

% data processing
par = pc(input);
par.prob_distro_function = 'Boltz';
par.tmesh_type = 'linear';
par = refresh_device(par);
xpos = 0;

for i = 1:length(var)
    par.B_ionic(4) = var(i);
    par = refresh_device(par);
    soleqs{i} = equilibrate(par);
end

for i = 1:length(soleqs)
    soleq = soleqs{i};
    sol = doCV(soleq.ion, 0, 0, -1, 1, 0.5, 1, 500);
    sols{i} = sol;
end

%% plot
df_plot.d2d_var(sols, xpos);
