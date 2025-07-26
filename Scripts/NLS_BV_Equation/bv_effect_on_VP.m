% - - - - - - - - - - CODE START - - - - - - - - - -

% pulse height and time
V_bias = 0; % stabilisation bias (0 V)
t_stab = 20; % stabilisation time

V_set = 0.5; % SET voltage
V_read = 0.1; % read voltage
V_pulse = [V_set, V_read, V_read, V_read, V_read, V_read]; % pulse voltage array
t_pulse_set = 5;
t_pulse_read = 5e-2;
t_pulse = [t_pulse_set, t_pulse_read, t_pulse_read, t_pulse_read, t_pulse_read, t_pulse_read];

% pulse shape
t_ramp = 1e-5; % ramp time for voltage transitions
t_cycle = 2; % one complete pulse time + relexation time
light_intensity = 0;

%% Simulation With BV

% parameter setting
% input = 'Input_files/single-layer.csv';
input = 'Input_files/mapi.csv';
par_bv = pc(input);
par_bv = refresh_device(par_bv);
soleq_bv = equilibrate(par_bv);

% transfer to simulation model
soleq_bv.eq.par.isEquilibrate = 'sim';
soleq_bv.ion.par.isEquilibrate = 'sim';

sol_bv = doPulse_v2(soleq_bv.ion, V_bias, V_pulse, t_pulse, t_cycle, t_stab, t_ramp, light_intensity);

%% Simulation Without BV

% parameter setting
par = pc(input);
par.j0 = 0;
par.mobseti = 0;
par = refresh_device(par);
soleq = equilibrate(par);

% transfer to simulation model
soleq.eq.par.isEquilibrate = 'sim';
soleq.ion.par.isEquilibrate = 'sim';

sol = doPulse_v2(soleq.ion, V_bias, V_pulse, t_pulse, t_cycle, t_stab, t_ramp, light_intensity);

%% Data Process
sols = cell(1, 2);
sols{1} = sol_bv;
sols{2} = sol;

%% Data Process
J_array = {};
t_array = {};

for i = 1:length(sols)

    t_offset = 0;

    for j = 1:(length(sols{i}) - 1)
        sol_temp = sols{i}{j + 1};

        [J_temp, ~, xmesh] = df_analysis.calcJ(sol_temp);
        xpos = 0;
        ppos = getpointpos(xpos, xmesh);
        J_peak = J_temp.tot(:, ppos);
        t_peak = sol_temp.t(:);

        % remove the outlier data points (by mean)
        mean_J = mean(J_peak); % mean value of J_peak
        % multiple = 5; % if the data is 5x larger than the mean_J => outlier
        % filter = (J_peak <= multiple * mean_J) & (J_peak >= 0);
        filter = J_peak >= 0;
        J_peak_filtered = J_peak(filter);
        t_peak_filtered = t_peak(filter);
        t_global = t_offset + t_peak_filtered; % global time series (for plot)

        % store the peak data for later sampling (reference point)
        J_array{i}{j} = J_peak_filtered;
        t_array{i}{j} = t_global;

        % t_global = t_offset + t_peak;
        % J_array{i}{j} = J_peak;
        % t_array{i}{j} = t_peak;

        % plot(t_global, J_peak_filtered, '.', 'LineWidth', 1);
        t_offset = t_offset + (2 * t_ramp) + t_pulse(j) + t_cycle;
    end

end
%%

figure('Name', 'Pulse Peak');
hold on;
set(gca, 'YScale', 'log');
colors = lines(2);

plot(NaN, NaN, '.', 'Color', colors(1,:), 'DisplayName', 'sol\_bv');
plot(NaN, NaN, '.', 'Color', colors(2,:), 'DisplayName', 'sol');

for i = 1:length(sols)

    for j = 1:length(J_array{i}) % Loop over each pulse in group
        plot(t_array{i}{j}, J_array{i}{j}, '.', 'Color', colors(i, :), 'HandleVisibility', 'off');
    end


end
legend('show');

%%
sampling_idx = 50;
J_sampling_pnts = {};
t_sampling_pnts = {};

for i = 1:length(J_array)

    J_sampling_pnts{i} = [];
    t_sampling_pnts{i} = [];

    for j = 2:length(J_array{i})

        if length(J_array{i}{j}) >= sampling_idx
            J_sampling_pnts{i}(end + 1) = J_array{i}{j}(sampling_idx);
            t_sampling_pnts{i}(end + 1) = t_array{i}{j}(sampling_idx);
        end

    end

end

%%
figure('Name', 'Sampling pnt'); 
hold on;
set(gca, 'YScale', 'log');

% for i = 1:length(J_sampling_pnts)
%     plot(t_sampling_pnts{i}, J_sampling_pnts{i}, '--o', 'LineWidth', 0.5);
% end

colors = lines(2);

for i = 1:length(J_sampling_pnts)

    if i == 1
        label = 'sol\_bv';
    else
        label = 'sol';
    end

    plot(t_sampling_pnts{i}, J_sampling_pnts{i}, ...
        '--o', 'LineWidth', 0.5, ...
        'Color', colors(i, :), ...
        'DisplayName', label);
end

legend('show');
