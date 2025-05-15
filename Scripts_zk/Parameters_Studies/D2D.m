
% Device-to-device (D2D) JV measurements,
% to compare the devices.

%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
inputs = {
    'Input_files_zk/peapi.csv';
    'Input_files_zk/peapi_if.csv';
    };

% - - - - - - - - - - handles
soleqs = cell(1, length(inputs));
sols = cell(1, length(inputs));
xpos = 0;
on_off_ratio_handle = cell(1, length(inputs));

%% - - - - - - - - - - DATA PROCESSING - - - - - - - - - -

for i = 1 : length(inputs)
    input = inputs{i};
    par = pc(input);
    par.prob_distro_function = 'Boltz';
    par.tmesh_type = 'linear';
    par = refresh_device(par);
    soleqs{i} = equilibrate(par);
end

%% - - - - - - - - - - DO MEASUREMENTS - - - - - - - - - -

for i = 1 : length(soleqs)
    soleq = soleqs{i};
    sol = doCV(soleq.ion, 0, 0, 5, -1, 5, 1, 100); % solution
    sols{i} = sol;
end

%% - - - - - - - - - - PLOTTING - - - - - - - - - -

figure('Name', 'Device-to-device Variability');
hold on;

for i = 1 : length(sols)
    xmesh = sols{i}.x;
    ppos = getpointpos(xpos, xmesh);

    J = dfana.calcJ(sols{i});
    Vapp = dfana.calcVapp(sols{i});
    t = sols{i}.t;
    
    % - - - - - - - - - - J - V Plot
    plot(Vapp, J.tot(:, ppos), 'DisplayName', ['D' num2str(i)]);
 
    % - - - - - - - - - - calculate and plot SET/RESET voltage point
    V_temp = Vapp;
    J_temp = J.tot(:, ppos);

    V_temp = V_temp(:); % reshape the data struc
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
    
    J_SET_fwd = J_temp(SET_fwd);
    J_SET_bwd = J_temp(SET_bwd);

    [~, J_read_out_idx_fwd] = min(abs(V_temp(SET_fwd) - V_read_out));
    J_read_out_fwd = J_SET_fwd(J_read_out_idx_fwd);
    
    [~, J_read_out_idx_bwd] = min(abs(V_temp(SET_bwd) - V_read_out));
    J_read_out_bwd = J_SET_bwd(J_read_out_idx_bwd);

    disp(['D2D.m: current values at V (fwd) = ', num2str(V_read_out), ' V: ', num2str(J_read_out_fwd)]);
    disp(['D2D.m: current values at V (bwd) = ', num2str(V_read_out), ' V: ', num2str(J_read_out_bwd)]);
   
    resistance_HRS = V_read_out / abs(J_read_out_bwd);
    resistance_LRS = V_read_out / abs(J_read_out_fwd);
    on_off_ratio = resistance_HRS / resistance_LRS;

    disp(['D2D.m: ON/OFF ratio: ', num2str(on_off_ratio)]);

    % - - - - - - - - - - save the ON/OFF ratio for each device
    on_off_ratio_handle{i}.device = i;
    on_off_ratio_handle{i}.ratio = on_off_ratio;
    
    % - - - - - - - - - - plot
    scatter(V_read_out, J_read_out_fwd, 50, 'rx', 'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    scatter(V_read_out, J_read_out_bwd, 50, 'rx', 'MarkerEdgeColor', 'k', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    plot([V_read_out, V_read_out], [J_read_out_fwd, J_read_out_bwd], 'k--', 'LineWidth', 0.75, 'HandleVisibility', 'off');
end
hold off

xlabel('Applied Voltage, Vapp [V]');
ylabel('Current Density, J [A cm^{-2}]');
legend('Location', 'northwest', 'FontSize', 10);

current_xlim = xlim; 
x_range = diff(current_xlim); 
increased_xlim = [current_xlim(1) - 0.125 * x_range, current_xlim(2) + 0.125 * x_range];
xlim(increased_xlim);

current_ylim = ylim;
y_range = diff(current_ylim);
increased_ylim = [current_ylim(1), current_ylim(1) + 1.25 * y_range];
ylim(increased_ylim);