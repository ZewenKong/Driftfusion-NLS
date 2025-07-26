% - - - - - - - - - - CODE START - - - - - - - - - -

%% Simulation

% parameter setting
input = 'Input_files/mapi.csv';
par = pc(input);
par.mu_c = [0, 0, 2.25e-12, 2.25e-12, 0];
par.mu_a = [0, 0, 0, 2.25e-12, 2.25e-12];
par.B = [6.3e-09, 0, 3.6e-10, 0, 6.8e-09];
par.j0 = 0;
par = refresh_device(par);

% do equilibrate
soleq = equilibrate(par);

% pulse height and time
V_bias = 0; % stabilisation bias (0 V)
t_stab = 10; % stabilisation time

V_set = 0.8; % SET voltage
V_read = 0.2; % read voltage
V_pulse = [V_set, V_read, V_read, V_read]; % pulse voltage array
t_pulse_set = 1e-1;
t_pulse_read = 1e-2;
t_pulse = [t_pulse_set, t_pulse_read, t_pulse_read, t_pulse_read, t_pulse_read];

% pulse shape
t_ramp = 1e-5; % ramp time for voltage transitions
t_cycle = 0.5; % one complete pulse time + relexation time
% tsample = [1e-1, 1e-2, 1e-2, 1e-2, 1e-2]; % time after pulse start to sample current
light_intensity = 0;

% transfer to simulation model
soleq.eq.par.isEquilibrate = 'sim';
soleq.ion.par.isEquilibrate = 'sim';

% do measurement
Pv2sol = doPulse_v2(soleq.ion, V_bias, V_pulse, t_pulse, t_cycle, t_stab, t_ramp, light_intensity);

%% Data Process & Plot
figure; hold on;
t_offset = 0;
J_array = {};
t_array = {};

for i = 1:(length(Pv2sol) - 1)
    sol_temp = Pv2sol{i + 1}; % create temp sol

    [J_temp, ~, xmesh] = df_analysis.calcJ(sol_temp);
    xpos = 0;
    ppos = getpointpos(xpos, xmesh);
    J_peak = J_temp.tot(:, ppos);
    t_peak = sol_temp.t(:);

    % remove the outlier data points (by mean)
    mean_J = mean(J_peak); % mean value of J_peak
    multiple = 2; % if the data is 5x larger than the mean_J => outlier
    filter = (J_peak <= multiple * mean_J) & (J_peak >= 0);

    J_peak_filtered = J_peak(filter);
    t_peak_filtered = t_peak(filter);
    t_global = t_offset + t_peak_filtered; % global time series (for plot)

    % store the peak data for later sampling (reference point)
    J_array{i} = J_peak_filtered;
    t_array{i} = t_global;

    plot(t_global, J_peak_filtered, '.', 'LineWidth', 1);
    t_offset = t_offset + (2 * t_ramp) + t_pulse(i) + t_cycle;
end

set(gca, 'YScale', 'log');

sampling_idx = 50;
J_sampling_pnts = [];
t_sampling_pnts = [];

for i = 1:length(J_array)

    if length(J_array{i}) >= sampling_idx
        J_sampling_pnts(end + 1) = J_array{i}(sampling_idx);
        t_sampling_pnts(end + 1) = t_array{i}(sampling_idx);
    end

end

plot(t_sampling_pnts, J_sampling_pnts, '--o', 'LineWidth', 0.5);
