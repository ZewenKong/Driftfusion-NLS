%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files/pmpi.csv';
active_thickness = [2.67e-5, 3.17e-5]; % var 1: active layer thickness
cation_mobility = [1e-9, 1e-10]; % var 2: cation mobility
var_1 = active_thickness;
var_2 = cation_mobility;

% - - - - - - - - - - handles
soleqs = cell(length(var_1), length(var_2));
sols = cell(size(soleqs));

% - - - - - - - - - - data processing
par = pc(input);
par.prob_distro_function = 'Boltz';
par.tmesh_type = 'linear';
par = refresh_device(par);
xpos = 0;

for i = 1:length(var_1)
    par.d(3) = var_1(i); % edit the value by accessing the position
    par = refresh_device(par); % refresh the par by driftfusion built-in func

    for j = 1:length(var_2)
        par.mu_c(3) = var_2(j);
        par = refresh_device(par);
        soleqs{i, j} = equilibrate(par);
    end

end

% - - - - - - - - - - do measurements
for i = 1:size(soleqs, 1)

    for j = 1:size(soleqs, 2)
        soleq = soleqs{i, j};
        sol = doCV(soleq.ion, 0, 0, 1, -1, 5e-1, 1, 500); % solution
        sols{i, j} = sol;
    end

end

%% - - - - - - - - - - plot
dfplot_ionic.d2d(sols, xpos, var_1, var_2);
