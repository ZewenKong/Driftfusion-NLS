function stats = CVstats(sol)
    % A function to pull statistics from a CV sweep using doCV
    % If multiple cycles have been performed, stats are taken from the first
    % cycle
    % sol - a solution from doCV

    % At present CVSTATS only works for P-I-N devices. Some adaption still
    % required for N-I-P devices

    % Check number of cycles
    num_cycles = sol.par.V_fun_arg(4);
    % Get the solution in a form of a sweep from V_min to V_max and back again.
    % This is complicated by the need to be able to accomodate people who don't
    % start their CV scan from 0 V
    if num_cycles > 1
        Vapp = df_analysis.calcVapp(sol);
        index = find(Vapp == sol.par.V_fun_arg(1), 2);
        one_sweep_index = index(2) - 1; % Index of last element of first cycle
        start = find(Vapp == min(Vapp), 1); % Index where V = V_min
        Vapp = Vapp(start:start + one_sweep_index)';
        sol.u = sol.u(start:start + one_sweep_index, :, :);
        sol.t = sol.t(start:start + one_sweep_index);
        J = df_analysis.calcJ(sol);
        J = J.tot(:, 1);
    else
        Vapp = df_analysis.calcVapp(sol);
        start = find(Vapp == min(Vapp), 1);
        J = df_analysis.calcJ(sol);
        J = J.tot(:, 1);
        J_half = J(start:end);
        J_other_half = J(2:start - 1);
        V_half = Vapp(start:end);
        V_other_half = Vapp(2:start - 1);

        if length(J_other_half) == 1
            J_half(end + 1) = J_other_half;
            J = J_half;
            V_half(end + 1) = V_other_half;
            Vapp = V_half;
        else
            J = [J_half; J_other_half];
            Vapp = [V_half V_other_half]';
        end

    end

    %% Define which data points correspond to forward and reverse scan directions
    change_sweep_direction_index = find(Vapp == max(Vapp), 1);
    J_f = J(1:change_sweep_direction_index);
    V_f = Vapp(1:change_sweep_direction_index);

    for i = 1:length(V_f)

        if abs(V_f(i)) < 1e-10
            V_f(i) = 0;
        end

    end

    J_r = J(change_sweep_direction_index:end);
    V_r = Vapp(change_sweep_direction_index:end);

    for i = 1:length(V_r)

        if abs(V_r(i)) < 1e-10
            V_r(i) = 0;
        end

    end

    %% Find stats for forward scan
    if length(J_f) == length(unique(J_f))
        stats.Jsc_f = interp1(V_f, J_f, 0, 'linear');
    else
        [Jtemp, rep_idx] = unique(J_f, 'stable', 'first');
        Vtemp = [V_f(1:rep_idx - 1), V_f(rep_idx + 1:end)];
        stats.Jsc_f = interp1(Vtemp, Jtemp, 0, 'linear');
    end

    if isnan(stats.Jsc_f)
        warning('No Jsc available- Vapp must pass through 0 to obtain Jsc')
        stats.Jsc_f = 0;
    end

    stats.Voc_f = interp1(J_f, V_f, 0, 'linear');

    if isnan(stats.Voc_f)
        warning('No Voc available- try increasing applied voltage range')
        stats.Voc_f = 0;
    end

    % Incident optical power
    Pin = df_analysis.calcPin(sol);

    A_f = 0; % Hysteresis Factor

    if stats.Jsc_f ~= 0 && stats.Voc_f ~= 0
        pow_f = J_f .* V_f;
        stats.mpp_f = min(pow_f);

        if stats.mpp_f > 0
            warning('Maximum power does not occur in the fourth quadrant')
            stats.mpp_f = 0;
            stats.mppV_f = 0;
        elseif stats.mpp_f < 0
            stats.mpp_f = abs(stats.mpp_f);
            stats.mppV_f = Vapp(-pow_f == stats.mpp_f);
        end

        stats.efficiency_f = 100 * (stats.mpp_f / Pin);
        stats.FF_f = -stats.mpp_f / (stats.Jsc_f * stats.Voc_f);

        try
            A_f = abs(trapz(V_f(V_f >= 0 & V_f <= stats.Voc_f), J_f(V_f >= 0 & V_f <= stats.Voc_f)));
        catch
            warning('Cannot calculate a hystersis factor, area may be too small.')
        end

    end

    %% Find stats for reverse scan
    stats.Jsc_r = interp1(V_r, J_r, 0, 'linear');

    if isnan(stats.Jsc_r)
        warning('No Jsc available- Vapp must pass through 0 to obtain Jsc')
        stats.Jsc_r = 0;
    end

    stats.Voc_r = interp1(J_r, V_r, 0, 'linear');

    if isnan(stats.Voc_r)
        warning('No Voc available- try increasing applied voltage range')
        stats.Voc_r = 0;
    end

    A_r = 0; %Hysteresis Factor

    if stats.Jsc_r ~= 0 && stats.Voc_r ~= 0
        pow_r = J_r .* V_r;
        stats.mpp_r = min(pow_r);

        if stats.mpp_r > 0
            warning('Maximum power does not occur in the fourth quadrant')
            stats.mpp_r = 0;
            stats.mppV_r = 0;
        elseif stats.mpp_r < 0
            stats.mpp_r = abs(stats.mpp_r);
            stats.mppV_r = Vapp(-pow_r == stats.mpp_r);
        end

        stats.efficiency_r = 100 * (stats.mpp_r / Pin);
        stats.mppV_r = Vapp(-pow_r == stats.mpp_r);
        stats.FF_r = -stats.mpp_r / (stats.Jsc_r * stats.Voc_r);

        try
            A_r = abs(trapz(V_r(V_r >= 0 & V_r <= stats.Voc_r), J_r(V_r >= 0 & V_r <= stats.Voc_r)));
        catch
            warning('Cannot calculate a hystersis factor, area may be too small.')
        end

    end

    %% Sign to identify inverted hysteresis
    if A_r ~= 0 && A_f ~= 0

        if A_r >= A_f
            B = 1;
        elseif A_r < A_f
            B = -1;
        end

        stats.HF = B * abs((A_r - A_f) / A_r);
    else
        stats.HF = 0;
    end

end
