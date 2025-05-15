
% Device-to-device (D2D) JV measurements,
% to compare the different properties (ion and vacancy mobility)

%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files_zk/peapi.csv';
ani_mu = [1e-12, 1e-11, 1e-10]; % iodine ions (in PCBM)
cat_mu = [1e-12, 1e-11, 1e-10]; % iodine vacancies (in RPP)
var_1 = ani_mu;
var_2 = cat_mu;

par = pc(input);
par.prob_distro_function = 'Boltz';
par.tmesh_type = 'linear';
par = refresh_device(par);
xpos = 0;

% - - - - - - - - - - handles
soleqs = cell(length(var_1), length(var_2));
sols = cell(size(soleqs));

%% - - - - - - - - - - DATA PROCESSING - - - - - - - - - -

for i = 1 : length(var_1)
    par.mu_a(4) = var_1(i); % anion mobility at interface
    par.mu_a(5) = var_1(i); % anion mobility in PCBM
    par = refresh_device(par); % refresh the device

    for j = 1 : length(var_2)
        par.mu_c(3) = var_2(i);
        par.mu_c(4) = var_2(i);
        par = refresh_device(par);

        soleqs{i, j} = equilibrate(par);
    end
end

%% - - - - - - - - - - DO MEASUREMENTS - - - - - - - - - -

for i = 1 : size(soleqs, 1)
    for j = 1 : size(soleqs, 2)

        soleq = soleqs{i, j};
        sol = doCV(soleq.ion, 0, 0, 3, -1, 5e-1, 1, 500);
        sols{i, j} = sol;
    end
end

%% - - - - - - - - - - PLOTTING - - - - - - - - - -

figure('Name', 'Anion/Cation Mobility Dependent JV');
hold on;
for i = 1 : size(sols, 1)
    for j = 1 : size(sols, 2)

        xmesh = sols{i, j}.x;
        ppos = getpointpos(xpos, xmesh);

        J = dfana.calcJ(sols{i, j});
        Vapp = dfana.calcVapp(sols{i, j});
        t = sols{i, j}.t;

        plot(Vapp, J.tot(:, ppos), 'DisplayName', [num2str(ani_mu(i)) '&' num2str(cat_mu(j))]);
    end
end
hold off;

xlabel('Applied Voltage, Vapp [V]');
ylabel('Current Density, J [A cm^{-2}]');
legend('show', 'Location', 'northwest', 'FontSize', 10);

current_xlim = xlim; 
x_range = diff(current_xlim); 
increased_xlim = [current_xlim(1) - 0.125 * x_range, current_xlim(2) + 0.125 * x_range];
xlim(increased_xlim);