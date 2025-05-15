
% Device-to-device (D2D) JV measurements,
% to compare the different properties (thickness and mobility).
%
% Define the ON/OFF ratio depencdence on two variables.

%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files_zk/mapi.csv';
active_thickness = [2.2e-5, 2.7e-5, 3.2e-5]; % var 1: active layer thickness
cation_mobility = [1e-8, 1e-10, 1e-12]; % var 2: cation mobility
var_1 = active_thickness;
var_2 = cation_mobility;

par = pc(input);
par.prob_distro_function = 'Boltz';
par.tmesh_type = 'linear';
par = refresh_device(par);
xpos = 0;

% - - - - - - - - - - handles
soleqs = cell(length(var_1), length(var_2));
sols = cell(size(soleqs));
on_off_ratio_handle = cell(size(sols));

%% - - - - - - - - - - DATA PROCESSING - - - - - - - - - -

for i = 1 : length(var_1)
    par.d(3) = var_1(i); % edit the value by accessing the position
    par = refresh_device(par); % refresh the par by driftfusion built-in func

    for j = 1 : length(var_2)
        par.mu_c(3) = var_2(j);
        par = refresh_device(par);
        soleqs{i, j} = equilibrate(par);
    end
end

%% - - - - - - - - - - DO MEASUREMENTS - - - - - - - - - -

for i = 1 : size(soleqs, 1)
    for j = 1 : size(soleqs, 2)
        soleq = soleqs{i, j};
        sol = doCV(soleq.ion, 0, 0, 1, -1, 5e-1, 1, 500); % solution
        sols{i, j} = sol;
    end
end

for i = 1 : size(sols, 1)
    for j = 1 : size(sols, 2)
        xmesh = sols{i, j}.x;
        ppos = getpointpos(xpos, xmesh);

        J = dfana.calcJ(sols{i, j});
        Vapp = dfana.calcVapp(sols{i, j});
        t = sols{i, j}.t;

        % - - - - - - - - - - calculate the ON/OFF ratio
        V_temp = Vapp;
        J_temp = J.tot(:, ppos);
        V_temp = V_temp(:);
        V_read_out = 1; % read-out voltage
        V_sweep_dirc = [0; diff(V_temp)];

        % - - - - - - - - - - define the SET forward/backward scan range
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
    
        disp(['D2D_props.m: current values at V (fwd) = ', num2str(V_read_out), ' V: ', num2str(J_read_out_fwd)]);
        disp(['D2D_props.m: current values at V (bwd) = ', num2str(V_read_out), ' V: ', num2str(J_read_out_bwd)]);

        % - - - - - - - - - - calculate out the HRS/LRS resistance
        resistance_HRS = V_read_out / abs(J_read_out_bwd);
        resistance_LRS = V_read_out / abs(J_read_out_fwd);

        % - - - - - - - - - - calculate out and save the ON/OFF ratio
        on_off_ratio = resistance_HRS / resistance_LRS;
        disp(['D2D_props.m: ON/OFF ratio: ', num2str(on_off_ratio)]);
        on_off_ratio_handle{i, j} = on_off_ratio;
    end
end

%% - - - - - - - - - - PLOTTING - - - - - - - - - -

on_off_matrix = cell2mat(on_off_ratio_handle);

[Xq, Yq] = meshgrid(linspace(min(var_1), max(var_1), 200), ...
    linspace(min(var_2), max(var_2), 200));

[X, Y] = meshgrid(var_1, var_2);

on_off_matrix_interp = griddata(X(:), Y(:), on_off_matrix(:), Xq, Yq, 'cubic');

figure('Name', 'ON/OFF Ratio Dependence Plot');
imagesc([min(var_1), max(var_1)], [min(var_2), max(var_2)], on_off_matrix_interp);
colorbar;
clim([min(on_off_matrix(:)), max(on_off_matrix(:))]);

xlabel('Active Layer Thickness (nm)');
ylabel('Cation Mobility (cm^2/Vs)');
xlim([min(var_1), max(var_1)]);
ylim([min(var_2), max(var_2)]);

colormap(turbo);
shading interp;