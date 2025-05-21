
% This script defines the scan rate (code from L.J.F.H).

%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files/peapi_v2.csv';
% scan_rts = logspace(-2, 2, 5);
scan_rts = [5, 6, 7, 8, 9, 10];

par = pc(input); 
par.prob_distro_function = 'Boltz'; 
par.tmesh_type = 'linear';
soleq = equilibrate(par);
xpos = 0;

% - - - - - - - - - - handle
sols = cell(1, length(scan_rts)); % solutions cell array

%% - - - - - - - - - - DATA PROCESSING - - - - - - - - - -

for i = 1 : length(scan_rts)
    sol = doCV(soleq.ion, 0, 0, -1, 1, scan_rts(i), 1, 500); % solution
    sols{i} = sol;
end

%% - - - - - - - - - - PLOTTING - - - - - - - - - -

figure('Name', 'Scan Rate Dependent JV');
hold on;
for i = 1 : length(sols)
    xmesh = sols{i}.x;
    ppos = getpointpos(xpos, xmesh);
    J = dfana.calcJ(sols{i});
    Vapp = dfana.calcVapp(sols{i});
    plot(Vapp, J.tot(:, ppos), 'DisplayName', [num2str(scan_rts(i)) ' Vs^{-1}'])
end
hold off;

xlabel('Applied Voltage, Vapp [V]');
ylabel('Current Density, J [A cm^{-2}]');
legend('show', 'Location', 'northwest', 'FontSize', 10);

current_xlim = xlim; 
x_range = diff(current_xlim); 
increased_xlim = [current_xlim(1) - 0.125 * x_range, current_xlim(2) + 0.125 * x_range];
xlim(increased_xlim);