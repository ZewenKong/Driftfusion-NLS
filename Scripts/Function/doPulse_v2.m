% - - - - - - - - - - CODE START - - - - - - - - - -

function Pv2sol = doPulse_v2(sol_ini, V_bias, V_pulse, tpulse, tcycle, tstab, tramp, light_intensity)

    % code from L.J.F.H
    % https://github.com/lucy-hart/Driftfusion (branch ESA)
    %
    % simulation of a stabilise and pulse (SaP) measurement

    % sol_ini, solution containing intitial conditions
    % V_bias, (array) of stabilisation biases
    % tstab, the time of the bias stabilised

    % V_pulse, array of voltages to sample at in the pulsed JV
    % tpulse, the length of the voltage pulse (time of a single pulse)
    % tcycle, the time between two pulses (the relaxation time)
    % tramp, the time to ramp (up/down)

    % light_intensity (0, dark condition)

    disp('doPulse_v2.m: pulse starting')

    num_bias = length(V_bias); % length of the V_bias
    num_pulses = length(V_pulse);

    % each row in Pv2sol is one V_bias;
    % 1st row is the stable result,
    % 2nd to nth are the V_pulse results
    Pv2sol = cell(num_bias, num_pulses + 1);

    for i = 1:num_bias % stabilise at V_bias

        if i == 1
            sol_ill = changeLight(sol_ini, light_intensity, 0, 1);
            par = sol_ill.par;
        else
            sol_ill = Pv2sol{i - 1, 1};
            par = Pv2sol{i - 1, 1}.par;
        end

        par.tmesh_type = 1;
        par.t0 = 0;
        par.tmax = 1e-2;
        par.tpoints = 100;
        par.V_fun_type = 'sweep';

        if i == 1
            par.V_fun_arg(1) = 0;
        else
            par.V_fun_arg(1) = V_bias(i - 1);
        end

        par.V_fun_arg(2) = V_bias(i);
        par.V_fun_arg(3) = 1e-2;

        sol = dfNLS(sol_ill, par); % call dfNLS.m

        par = sol.par;
        par.tmesh_type = 1;
        par.t0 = 0;
        par.tmax = tstab; % hold device at V_bias for t_stab
        par.tpoints = 100;
        par.V_fun_type = 'constant';
        par.V_fun_arg(1) = V_bias(i);

        disp(['doPulse_v2.m: stabilising solution at ' num2str(V_bias(i)) ' V']);
        Pv2sol{i, 1} = dfNLS(sol, par);
    end

    for i = 1:num_bias % perform the pulsed JV for each V_bias (stabilise at V_bias)
        disp(['doPulse_v2.m: starting pulse for V_stab = ' num2str(V_bias(i)) ' V']);

        for j = 1:num_pulses % apply pulse voltage to every stabilised bias

            if j == 1
                sol = Pv2sol{i, 1}; % if j = 1, ramping up from background
                par = Pv2sol{i, 1}.par;
            else
                sol = sol; % if j \= 1, inheritance
                par = sol.par;
            end

            % - - - - - - - - - - VOLTAGE RAMPING UP - - - - - - - - - -

            par.mobseti = 1; % ion motion switch
            par.tmesh_type = 1;
            par.t0 = 0;
            par.tmax = tramp; % ramping up time taken
            par.tpoints = 100;

            par.V_fun_type = 'sweep';
            par.V_fun_arg(1) = V_bias(i);
            par.V_fun_arg(2) = V_pulse(j);
            par.V_fun_arg(3) = tramp;

            try
                sol = dfNLS(sol, par);
            catch
                warning(['doPulse_v2.m: could not ramp to voltage for V_pulse = ' num2str(V_pulse(j)) ' V'])
                sol = 0;
            end

            if not(isa(sol, 'struct'))
                Pv2sol{i, j + 1}.J_pulse = 0;
            elseif isa(sol, 'struct')

                % - - - - - - - - - - VOLTAGE PLATEAU - - - - - - - - - -

                par = sol.par;
                par.mobseti = 1; % ion motion switch
                par.tmesh_type = 1;

                par.t0 = 0;
                par.tmax = tpulse(j); % tpulse, length of the pulse peak
                par.tpoints = 100;

                par.V_fun_type = 'constant';
                par.V_fun_arg(1) = V_pulse(j);

                disp(['doPulse_v2.m: V_pulse = ' num2str(V_pulse(j)) ' V']);

                sol = dfNLS(sol, par); % assign the solution
                Pv2sol{i, j + 1} = sol;

                % [J_peak, ~, xmesh] = df_analysis.calcJ(Pv2sol{i, j + 1});
                % xpos = 0;
                % ppos = getpointpos(xpos, xmesh);
                % J_peak = J_peak.tot(:, ppos);
                % Pv2sol{i, j + 1}.J_peak = J_peak; % total current during pulse peak
                % Pv2sol{i, j + 1}.t_peak = Pv2sol{i, j + 1}.t; % pulse peak time

                % - - - - - - - - - - VOLTAGE RAMPING DOWN - - - - - - - - - -

                par.tmesh_type = 1;
                par.t0 = 0;
                par.tmax = tramp;
                par.tpoints = 100;

                par.V_fun_type = 'sweep';
                par.V_fun_arg(1) = V_pulse(j);
                par.V_fun_arg(2) = V_bias(i);
                par.V_fun_arg(3) = tramp;

                try
                    sol = dfNLS(Pv2sol{i, j + 1}, par);
                catch
                    warning(['doPulse_v2.m: could not ramp down to V_bias for V_pulse = ' num2str(V_pulse(j)) ' V'])
                    sol = Pv2sol{i, j + 1};
                end

                % - - - - - - - - - - RELAXATION - - - - - - - - - -

                par.tmesh_type = 1;
                par.t0 = 0;
                par.tmax = tcycle;
                par.tpoints = 100;

                par.V_fun_type = 'constant';
                par.V_fun_arg(1) = V_bias(i);

                try
                    sol = dfNLS(sol, par);
                catch
                    warning('doPulse_v2.m: could not go back to V_bias')
                    sol = Pv2sol{i, j + 1};
                end

            end

        end

    end

end
