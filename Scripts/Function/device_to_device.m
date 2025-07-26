% - - - - - - - - - - CODE START - - - - - - - - - -

%% Simulation

% parameter setting
inputs = {
          'Input_files/pmpi.csv';
          'Input_files/mapi.csv';
          };

% handles
soleqs = cell(1, length(inputs));
sols = cell(1, length(inputs));

for i = 1:length(inputs)
    % parameter setting
    input = inputs{i};
    par = pc(input);
    par.prob_distro_function = 'Boltz';
    par = refresh_device(par);

    % do equilibrium
    soleqs{i} = equilibrate(par);
end

% do measurements
for i = 1:length(soleqs)
    soleq = soleqs{i};
    sol = doCV(soleq.ion, 0, 0, -1, 1, 1e-1, 1, 500); % solution
    sols{i} = sol;
end

%% plot

% plot setting
xpos = 0;

% do plot
df_plot.d2d_var(sols, xpos);
