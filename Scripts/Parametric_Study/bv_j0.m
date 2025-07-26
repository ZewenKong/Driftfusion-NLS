%% - - - - - - - - - - CODE START - - - - - - - - - -

% data inputs
input = 'Input_files/pmpi_v2.csv';
j0s = [1e-3, 1e-2, 1e-1, 0, 1, 10, 100, 1000];
var = j0s;

% handles
soleqs = cell(1, length(var));
sols = cell(1, length(soleqs)); % solutions cell array

% data processing
par = pc(input);
par.prob_distro_function = 'Boltz';
par.tmesh_type = 'linear';
par = refresh_device(par);
xpos = 0;

for i = 1:length(var)
    par.j0 = var(i);
    par = refresh_device(par);
    soleqs{i} = equilibrate(par);
    soleqs{i}.el.par.isEquilibrate = "sim";
    soleqs{i}.ion.par.isEquilibrate = "sim";
end

%% do measurements
for i = 1:length(soleqs)
    soleq = soleqs{i};
    sol = doCV(soleq.ion, 0, 0, -1, 1, 0.05, 1, 500); % solution
    sols{i} = sol;
end

%% plot
df_plot.d2d_var(sols, xpos);
