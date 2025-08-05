function Pv1sol = doPulse_v1(sol_ini, V_bias, V_pulse, tpulse, tcycle, tstab, light_intensity)

    disp('doPulse_v1.m: starting SaP')

    num_bias = length(V_bias);
    num_pulses = length(V_pulse);

    Pv1sol = cell(num_bias, num_pulses + 1);
    sol_ill = changeLight(sol_ini, light_intensity, 0, 1);

    for i = 1:num_bias

        par = sol_ill.par;
        par.tmesh_type = 1;
        par.t0 = 0;
        par.tmax = 1e-2;
        par.tpoints = 100;
        par.V_fun_type = 'sweep';
        par.V_fun_arg(1) = 0;
        par.V_fun_arg(2) = V_bias(i);
        par.V_fun_arg(3) = 1e-2;

        sol = dfII(sol_ill, par);
        par = sol.par;

        par.tmesh_type = 1;
        par.t0 = 0;
        par.tmax = tstab;
        par.tpoints = 100;
        par.V_fun_type = 'constant';
        par.V_fun_arg(1) = V_bias(i);

        disp(['doPulse_v1.m: stabilising solution at ' num2str(V_bias(i)) ' V'])
        Pv1sol{i, 1} = dfII(sol, par);

    end

    for i = 1:num_bias
        disp(['doPulse_v1.m: starting doPulse for V_stab = ' num2str(V_bias(i)) ' V'])
        par = Pv1sol{i, 1}.par;
        par.mobseti = 1;

        par.tmesh_type = 1;
        par.t0 = 0;
        par.tmax = tcycle;

        % par.tpoints = 1e2;
        max_dt = 1e-4;
        par.tpoints = ceil((par.tmax - par.t0) / max_dt) + 1;

        par.V_fun_type = 'smoothed_square';
        par.V_fun_arg(1) = V_bias(i);

        for j = 1:num_pulses

            par.V_fun_arg(2) = V_pulse(j);
            par.V_fun_arg(3) = tcycle;
            duty_cycle = 100 * tpulse(j) / tcycle;
            par.V_fun_arg(4) = duty_cycle; % the ratio of pulse time in total time
            disp(['doPulse_v1.m: V_pulse = ' num2str(V_pulse(j)) ' V'])
            Pv1sol{i, j + 1} = dfII(Pv1sol{i, 1}, par);

        end

    end

end
