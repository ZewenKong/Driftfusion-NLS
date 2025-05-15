
% Device-to-device (D2D) JV measurements,
% to compare the different property (thickness here).

%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files_zk/peapi.csv';
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

for i = 1 : length(var)
    par.d(3) = var(i); % active layer thickness
    par = refresh_device(par);
    soleqs{i} = equilibrate(par);
end

%% - - - - - - - - - - DO MEASUREMENTS - - - - - - - - - -

for i = 1 : length(soleqs)
    soleq = soleqs{i};
    sol = doCV(soleq.ion, 0, 0, 1, -1, 5e-1, 1, 500); % solution
    sols{i} = sol;
end

%% - - - - - - - - - - PLOTTING - - - - - - - - - -

figure('Name', 'Device-to-device Variability on Physical Property');
hold on;
for i = 1 : length(sols)
    xmesh = sols{i}.x;
    ppos = getpointpos(xpos, xmesh);

    J = dfana.calcJ(sols{i});
    Vapp = dfana.calcVapp(sols{i});
    
    % - - - - - - - - - - J - V Plot
    plot(Vapp, J.tot(:, ppos), 'DisplayName', ['t' num2str(var(i))]);

    % - - - - - - - - - - calculate and plot SET/RESET voltage point
    V_temp = Vapp;
    J_temp = J.tot(:, ppos);

    V_temp = V_temp(:); % reshape the data struc
    J_temp = J_temp(:);
    
    dJdV = gradient(J_temp) ./ gradient(V_temp); % gradient of each point
    
    % - - - - - - - - - - define SET/RESET range
    V_sweep_dirc = [0; diff(V_temp)]; % check the scan direction
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
    scatter(V_SET, J_SET, 50, 'rx', 'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    scatter(V_RESET, J_RESET, 50, 'bx', 'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    
    % - - - - - - - - - - calculate the ON/OFF ratio
    V_temp = Vapp;
    J_temp = J.tot(:, ppos);
    V_temp = V_temp(:);

    V_read_out = 1;
    V_sweep_dirc = [0; diff(V_temp)];

    % - - - - - - - - - - define the forward/backward scan
    SET_fwd = (V_temp > 0) & (V_sweep_dirc > 0);
    SET_bwd = (V_temp > 0) & (V_sweep_dirc < 0);
    
    % - - - - - - - - - - find the current value in the SET forward/backward scan range
    J_SET_fwd = J_temp(SET_fwd);
    J_SET_bwd = J_temp(SET_bwd);

    % - - - - - - - - - - find the most closest voltage point (close to read out voltage)
    [~, J_read_out_idx_fwd] = min(abs(V_temp(SET_fwd) - V_read_out));
    J_read_out_fwd = J_SET_fwd(J_read_out_idx_fwd);
    
    [~, J_read_out_idx_bwd] = min(abs(V_temp(SET_bwd) - V_read_out));
    J_read_out_bwd = J_SET_bwd(J_read_out_idx_bwd);

    disp(['D2D_prop.m: current values at V (fwd) = ', num2str(V_read_out), ' V: ', num2str(J_read_out_fwd)]);
    disp(['D2D_prop.m: current values at V (bwd) = ', num2str(V_read_out), ' V: ', num2str(J_read_out_bwd)]);
    
    % - - - - - - - - - - calculate out the HRS/LRS resistance
    resistance_HRS = V_read_out / abs(J_read_out_bwd);
    resistance_LRS = V_read_out / abs(J_read_out_fwd);

    % - - - - - - - - - - calculate out and save the ON/OFF ratio
    on_off_ratio = resistance_HRS / resistance_LRS;
    disp(['D2D_prop.m: ON/OFF ratio: ', num2str(on_off_ratio)]);

    on_off_ratio_handle{i}.device = i;
    on_off_ratio_handle{i}.ratio = on_off_ratio;
    
    % - - - - - - - - - - visualise the the ON/OFF ratio
    scatter(V_read_out, J_read_out_fwd, 50, 'rx', 'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    scatter(V_read_out, J_read_out_bwd, 50, 'rx', 'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    plot([V_read_out, V_read_out], [J_read_out_fwd, J_read_out_bwd], 'k--', 'LineWidth', 0.75, 'HandleVisibility', 'off');
end
hold off;

xlabel('Applied Voltage, Vapp [V]');
ylabel('Current Density, J [A cm^{-2}]');
legend('Location', 'northwest', 'FontSize', 10);

current_ylim = ylim;
y_range = diff(current_ylim);
increased_ylim = [current_ylim(1), current_ylim(1) + 1.25 * y_range];
ylim(increased_ylim);