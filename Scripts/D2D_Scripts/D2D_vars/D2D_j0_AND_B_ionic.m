%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files/pmpi.csv';

j0_array = [1e-12, 5e-12]; % j0 numerical array, double row vector (dbv)
B_ionic_array = [1e-15, 5e-15]; % B_ionic double row vector
var_1 = j0_array;
var_2 = B_ionic_array;
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
    par.j0 = var_1(i); % j0
    par = refresh_device(par);

    for j = 1:length(var_2)
        par.B_ionic(4) = var_2(j);
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
dfplot_ionic.d2d(sols, xpos, var_1, var_2);
