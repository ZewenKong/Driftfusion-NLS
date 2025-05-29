% This script defines the scan rate (code from L.J.F.H).

%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files/pmpi_v2.csv';
scan_rts = [1e-1, 5e-1, 1, 5];

par = pc(input);
par.prob_distro_function = 'Boltz';
par.tmesh_type = 'linear';
soleq = equilibrate(par);

% - - - - - - - - - - handle
sols = cell(1, length(scan_rts)); % solutions cell array

%% - - - - - - - - - - DATA PROCESSING - - - - - - - - - -

for i = 1:length(scan_rts)
    sol = doCV(soleq.ion, 0, 0, -1, 1, scan_rts(i), 1, 500); % solution
    sols{i} = sol;
end

%% - - - - - - - - - - PLOTTING - - - - - - - - - -

figure('Name', 'Scan Rate Dependent JV');
hold on;

for i = 1:length(sols)
    xmesh = sols{i}.x;
    xpos = 0; % position
    ppos = getpointpos(xpos, xmesh);
    J = dfana.calcJ(sols{i});
    Vapp = dfana.calcVapp(sols{i});
    plot(Vapp, J.tot(:, ppos), 'DisplayName', [num2str(scan_rts(i)) ' Vs^{-1}'])
end

hold off;

xlabel('Applied Voltage, Vapp [V]');
ylabel('Current Density, J [A cm^{-2}]');
legend('show', 'Location', 'northwest', 'FontSize', 10);
