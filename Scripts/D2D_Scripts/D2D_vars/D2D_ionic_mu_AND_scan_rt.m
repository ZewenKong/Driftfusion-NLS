%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files/pmpi_v2.csv';

ionic_mus = [1e-12, 1e-14, 1e-16, 1e-18];
scan_rts = [5e-2, 7.5e-2, 1e-1, 1.25e-1];
var_1 = ionic_mus;
var_2 = scan_rts;
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
    par.mu_c(3) = var_1(i);
    par.mu_c(4) = var_1(i);
    par.mu_a(4) = var_1(i);
    par.mu_a(5) = var_1(i);
    par = refresh_device(par);
    soleqs{i} = equilibrate(par);
end

% - - - - - - - - - - do measurements
for i = 1:length(soleqs)
    soleq = soleqs{i};

    for j = 1:length(var_2)
        sol = doCV(soleq.ion, 0, 0, -1, 1, var_2(j), 1, 500); % solution
        sols{j} = sol;
    end

end

%% - - - - - - - - - - plot
dfplot_ionic.d2d_var(sols, xpos);
