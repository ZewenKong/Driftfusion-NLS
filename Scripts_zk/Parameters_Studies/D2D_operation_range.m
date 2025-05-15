
% Device-to-device (D2D) JV measurements,
% defines the range of bias voltages to apply in order to 
% produce a significant change in current density.

%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files_zk/peapi.csv';
par = pc(input);
par.prob_distro_function = 'Boltz';
par.tmesh_type = 'linear';
par = refresh_device(par);

cycle = 1; % cycle value
maxV = [5, 10]; % max V value
minV = 0; % min V value (consider the SET process)
xpos = 0;

% - - - - - - - - - - handle
sols = cell(1, length(maxV)); % solutions cell array

%% - - - - - - - - - - DATA PROCESSING - - - - - - - - - -

soleq = equilibrate(par);

%% - - - - - - - - - - DO MEASUREMENTS - - - - - - - - - -

for i = 1 : length(maxV)
    sol = doCV(soleq.ion, 0, 0, maxV(i), minV, 5, cycle, 100); % solution
    sols{i} = sol; % add the solution to the cell array
end

%% - - - - - - - - - - PLOTTING - - - - - - - - - -

figure('Name', 'Operation Range C.F.');
hold on;
for i = 1 : length(sols)
    xmesh = sols{i}.x;
    ppos = getpointpos(xpos, xmesh);
    J = dfana.calcJ(sols{i}); % calculate the J of the solution
    Vapp = dfana.calcVapp(sols{i}); % calculate the V of the solution
    t = sols{i}.t; % time
    plot(t, J.tot(:, ppos), 'DisplayName', [num2str(maxV(i)) 'V']);
end
hold off;

xlabel('Time, t [s]');
ylabel('Current Density, J [A cm^{-2}]');
legend('Location', 'northwest', 'FontSize', 10);

current_xlim = xlim; 
x_range = diff(current_xlim); 
increased_xlim = [0.1, current_xlim(2) + 0.125 * x_range];
xlim(increased_xlim);

current_ylim = ylim;
y_range = diff(current_ylim);
increased_ylim = [current_ylim(1), current_ylim(1) + 1.25 * y_range];
ylim(increased_ylim);