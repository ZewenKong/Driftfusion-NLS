% Device-to-device (D2D) JV measurements

%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
inputs = {
          'Input_files/pmpi_v2.csv';
          'Input_files/pmpi_uni.csv';
          };
xpos = 0;

% - - - - - - - - - - handles
soleqs = cell(1, length(inputs));
sols = cell(1, length(inputs));

% - - - - - - - - - - data processings
for i = 1:length(inputs)
    input = inputs{i};
    par = pc(input);
    par.prob_distro_function = 'Boltz';
    par.tmesh_type = 'linear';
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
