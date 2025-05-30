%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files/pmpi_v2.csv';
ani_mu = [1e-9, 1e-10, 1e-11]; % iodine ions (in PCBM)
cat_mu = [1e-9, 1e-10, 1e-11]; % iodine vacancies (in RPP)
var_1 = ani_mu;
var_2 = cat_mu;
xpos = 0;

% - - - - - - - - - - handles
soleqs = cell(length(var_1), length(var_2));
sols = cell(size(soleqs));

% - - - - - - - - - - data processing
par = pc(input);
par.prob_distro_function = 'Boltz';
par.tmesh_type = 'linear';
par = refresh_device(par);

for i = 1:length(var_1)
    par.mu_a(4) = var_1(i); % anion mobility at interface
    par.mu_a(5) = var_1(i); % anion mobility in PCBM
    par = refresh_device(par); % refresh the device

    for j = 1:length(var_2)
        par.mu_c(3) = var_2(i);
        par.mu_c(4) = var_2(i);
        par = refresh_device(par);
        soleqs{i, j} = equilibrate(par);
    end

end

% - - - - - - - - - - do measurements
for i = 1:size(soleqs, 1)

    for j = 1:size(soleqs, 2)
        soleq = soleqs{i, j};
        sol = doCV(soleq.ion, 0, 0, -1, 1, 1e-1, 1, 500);
        sols{i, j} = sol;
    end

end

%% - - - - - - - - - - plot
dfplot_ionic.d2d_vars(sols, xpos, var_1, var_2);
