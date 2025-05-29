% Device-to-device (D2D) JV measurements,
% to compare the devices.

%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
inputs = {
          'Input_files/pmpi_v2.csv';
          'Input_files/pmpi_uni.csv';
          };

% - - - - - - - - - - handles
soleqs = cell(1, length(inputs));
sols = cell(1, length(inputs));

%% - - - - - - - - - - DATA PROCESSING - - - - - - - - - -

for i = 1:length(inputs)
    input = inputs{i};
    par = pc(input);
    par.prob_distro_function = 'Boltz';
    par.tmesh_type = 'linear';
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

figure('Name', 'd2d');
hold on;

for i = 1:length(sols)
    xmesh = sols{i}.x;
    xpos = 0;
    ppos = getpointpos(xpos, xmesh);

    J = dfana.calcJ(sols{i});
    Vapp = dfana.calcVapp(sols{i});
    t = sols{i}.t;

    % - - - - - - - - - - J - V Plot
    plot(Vapp, J.tot(:, end), 'DisplayName', ['D' num2str(i)]);
end

hold off

xlabel('Applied Voltage, Vapp [V]');
ylabel('Current Density, J [A cm^{-2}]');
legend('Location', 'northwest', 'FontSize', 10);
