
% Cycle-to-cycle (C2C) JV measurement,
% to define the endurance and the stochasticity

% - - - - - - - - - - initialisation
equilibrate_init();
initialise_df;

%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - read in data files
input = 'Input_files_zk/peapi.csv';

par = pc(input); % original parameters
par.tmesh_type = 'linear'; % log10

par = refresh_device(par); % refresh the device

soleq = equilibrate(par);
cycle = 1; % cycle value

%% - - - - - - - - - - DO MEASUREMENT - - - - - - - - - -

sol = doCV(soleq.ion, 0, 0, -1, 1, 5e-2, cycle, 500); % solution

xmesh = sol.x; xpos = 0;
ppos = getpointpos(xpos, xmesh);

J = dfana.calcJ(sol);
Vapp = dfana.calcVapp(sol);
t = sol.t;

%% - - - - - - - - - - PLOTTING - - - - - - - - - -

figure('Name', 'Cycle-to-cycle Variability');
hold on;

% assume each cycle has same number of data points get the number of data points per cycle
data_pnts_per_cyc = length(Vapp) / (cycle);

%  - - - - - - - - - - handles
legend_handle = [];
legend_label = cell(1, cycle);
on_off_ratio_handle = cell(1, cycle); % ON/OFF ratio cell array (handle)

for i = 1 : cycle

    % round (), make sure idx is an integer
    % (1 : data_pnts_per_cyc), current data points array
    % e.g., round((1 - 1) * (data_pnts_per_cyc), 1st Cycle starts from 0
    %       round((2 - 1) * (data_pnts_per_cyc), 2nd Cycle starts after 1st Cycle

    cycle_idx = round(round(i - 1) * (data_pnts_per_cyc) + (1 : data_pnts_per_cyc));

    % - - - - - - - - - - J - V Plot
    plots = plot(Vapp(cycle_idx), J.tot(cycle_idx, ppos), 'LineWidth', 0.75);

    legend_handle = [legend_handle, plots];
    legend_label{i} = ['C' num2str(i)];

    % - - - - - - - - - - calculate and plot SET/RESET voltage point
    V_temp = Vapp(cycle_idx);
    J_temp = J.tot(cycle_idx, ppos);

    V_temp = V_temp(:); % re-shape the data struc
    J_temp = J_temp(:);

    dJdV = gradient(J_temp) ./ gradient(V_temp); % gradient of each point

    % - - - - - - - - - - define SET/RESET range
    V_sweep_dirc = [0; diff(V_temp)];
    SET_V_range = (V_temp > 0) & (V_sweep_dirc > 0);
    RESET_V_range = (V_temp < 0) & (V_sweep_dirc < 0); 

    % - - - - - - - - - - find the max gradient in SET range
    [max_slope_in_SET, idx_in_SET] = max(dJdV(SET_V_range));
    V_temp_in_SET = V_temp(SET_V_range);
    V_SET = V_temp_in_SET(idx_in_SET);
    J_temp_in_SET = J_temp(SET_V_range);
    J_SET = J_temp_in_SET(idx_in_SET);

    % - - - - - - - - - - find the max gradient in RESET range
    [max_slope_in_RESET, idx_in_RESET] = max(dJdV(RESET_V_range));
    V_temp_in_RESET = V_temp(RESET_V_range);
    V_RESET = V_temp_in_RESET(idx_in_RESET);
    J_temp_in_RESET = J_temp(RESET_V_range);
    J_RESET = J_temp_in_RESET(idx_in_RESET);

    % - - - - - - - - - - plot the SET/RESET points for each cycle
    scatter(V_SET, J_SET, 50, 'rx', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    scatter(V_RESET, J_RESET, 50, 'bx', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

    % - - - - - - - - - - calculate the ON/OFF ratio
    V_temp = Vapp(cycle_idx);
    J_temp = J.tot(cycle_idx, ppos);
    V_temp = V_temp(:);

    V_read_out = 1; % read-out voltage
    V_sweep_dirc = [0; diff(V_temp)];

    % - - - - - - - - - - define the forward/backward SET scan
    SET_fwd = (V_temp > 0) & (V_sweep_dirc > 0);
    SET_bwd = (V_temp > 0) & (V_sweep_dirc < 0);

    J_SET_fwd = J_temp(SET_fwd);
    J_SET_bwd = J_temp(SET_bwd);

    [~, J_read_out_idx_fwd] = min(abs(V_temp(SET_fwd) - V_read_out));
    J_read_out_fwd = J_SET_fwd(J_read_out_idx_fwd);

    [~, J_read_out_idx_bwd] = min(abs(V_temp(SET_bwd) - V_read_out));
    J_read_out_bwd = J_SET_bwd(J_read_out_idx_bwd);

    disp(['C2C.m: current values at V (fwd) = ', num2str(V_read_out), ' V: ', num2str(J_read_out_fwd)]);
    disp(['C2C.m: current values at V (bwd) = ', num2str(V_read_out), ' V: ', num2str(J_read_out_bwd)]);

    resistance_HRS = V_read_out / abs(J_read_out_bwd);
    resistance_LRS = V_read_out / abs(J_read_out_fwd);
    on_off_ratio = resistance_HRS / resistance_LRS;

    % - - - - - - - - - - save the ON/OFF ratio for each cycle
    on_off_ratio_handle{i}.cycle = i;
    on_off_ratio_handle{i}.ratio = on_off_ratio;
    disp(['C2C.m: ON/OFF ratio: ', num2str(on_off_ratio)]);

    % - - - - - - - - - - plot
    scatter(V_read_out, J_read_out_fwd, 50, 'rx', 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
    scatter(V_read_out, J_read_out_bwd, 50, 'rx', 'MarkerEdgeColor', 'k', 'LineWidth', 0.75);
    plot([V_read_out, V_read_out], [J_read_out_fwd, J_read_out_bwd], '--', 'LineWidth', 0.75);
end
hold off;
box on;

xlabel('Applied Voltage, Vapp [V]');
ylabel('Current Density, J [A cm^{-2}]');
legend(legend_handle, legend_label, 'Location', 'northwest', 'FontSize', 10);

% - - - - - - - - - - adjust the X-axis limit
current_xlim = xlim;
x_range = diff(current_xlim);
increased_xlim = [current_xlim(1) - 0.125 * x_range, current_xlim(2) + 0.125 * x_range];
xlim(increased_xlim);

% - - - - - - - - - - adjust the Y-axis limit
current_ylim = ylim;
y_range = diff(current_ylim);
increased_ylim = [current_ylim(1), current_ylim(1) + 1.25 * y_range];
ylim(increased_ylim);