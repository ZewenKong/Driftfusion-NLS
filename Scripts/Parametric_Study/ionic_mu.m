%% - - - - - - - - - - CODE START - - - - - - - - - -

input = 'Input_files/pmpi_v2.csv';
ionic_mus = [1e-10, 1e-11, 1e-12, 1e-13, 1e-14];
var = ionic_mus;

% handles
soleqs = cell(1, length(var));
sols = cell(size(soleqs));

% data processing
par = pc(input);
par.prob_distro_function = 'Boltz';
par.tmesh_type = 'linear';
par = refresh_device(par);
xpos = 0;

for i = 1:length(var)
    par.mu_c(3) = var(i);
    par.mu_c(4) = var(i);
    par.mu_a(4) = var(i);
    par.mu_a(5) = var(i);
    par = refresh_device(par);
    soleqs{i} = equilibrate(par);
    soleqs{i}.el.par.isEquilibrate = "sim";
    soleqs{i}.ion.par.isEquilibrate = "sim";
end

% do measurements
for i = 1:length(soleqs)
    soleq = soleqs{i};
    sol = doCV(soleq.ion, 0, 0, -1, 1, 1e-1, 1, 500); % solution
    sols{i} = sol;
end

%% plot
df_plot.d2d_var(sols, xpos);
