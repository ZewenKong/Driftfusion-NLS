% Device-to-device (D2D) JV measurements,
% to compare the different property (thickness here).

%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files/pmpi.csv';
active_thickness = [2.2e-5, 2.7e-5, 3.2e-5];
var = active_thickness;

par = pc(input);
par.prob_distro_function = 'Boltz';
par.tmesh_type = 'linear';
par = refresh_device(par);
xpos = 0;

% - - - - - - - - - - handles
soleqs = cell(1, length(var));
sols = cell(size(soleqs));
on_off_ratio_handle = cell(size(sols));

%% - - - - - - - - - - DATA PROCESSING - - - - - - - - - -

for i = 1:length(var)
    par.d(3) = var(i); % active layer thickness
    par = refresh_device(par);
    soleqs{i} = equilibrate(par);
end

%% - - - - - - - - - - DO MEASUREMENTS - - - - - - - - - -

for i = 1:length(soleqs)
    soleq = soleqs{i};
    sol = doCV(soleq.ion, 0, 0, -1, 1, 1e-1, 1, 500); % solution
    sols{i} = sol;
end

%% - - - - - - - - - - PLOTTING - - - - - - - - - -

figure('Name', 'Device-to-device Variability on Physical Property');
hold on;

for i = 1:length(sols)
    xmesh = sols{i}.x;
    xpos = 0;
    ppos = getpointpos(xpos, xmesh);

    J = dfana.calcJ(sols{i});
    Vapp = dfana.calcVapp(sols{i});

    plot(Vapp, J.tot(:, ppos), 'DisplayName', ['thickness = ' num2str(var(i))]);
end

hold off;

xlabel('Applied Voltage, Vapp [V]');
ylabel('Current Density, J [A cm^{-2}]');
legend('Location', 'northwest', 'FontSize', 10);
