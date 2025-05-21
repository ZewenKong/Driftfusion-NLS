function sol_CV = doCV(sol_ini, light_intensity, V0, Vmax, Vmin, scan_rate, cycles, tpoints)

    % doCV(ion solution, light intensity, starting voltage, max voltage,
    %      min voltage, scan rate, cycles, time points);
    %
    % Performs a cyclic voltammogram (CV) simulation
    % Input arguments:
    % SOL_INI = solution containing intitial conditions
    % LIGHT_INTENSITY = Light intensity for bias light (Suns)
    % V0 = Starting voltage (V)
    % VMAX = Maximum voltage point (V)
    % VMIN = Minimum voltage point (V)
    % SCAN_RATE = Scan rate (Vs-1)
    % CYCLES = No. of scan cycles
    % TPOINTS = No. of points in output time array
    %
    % LICENSE
    % Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
    % Imperial College London
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU Affero General Public License as published
    % by the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    %
    %% - - - - - - - - - - CODE START - - - - - - - - - -

    disp('doCV.m: starting CV scan')

    par = sol_ini.par;

    if light_intensity > 0 % set light intensity
        sol = lightonRs(sol_ini, light_intensity, -1, 0, par.Rs, 10);
    else
        sol = sol_ini;
    end

    % switch to V0 if different from initial conditions,
    % trying to do this within the triangular wave function resulted in
    % convergence issues so safer to do this here

    if V0 ~= par.V_fun_arg(1)
        sol = genVappStructs(sol, V0, 0);
    end

    % calculate tmax from scan rate and absolute change in voltage, deltaV

    deltaV = abs(Vmax - V0) + abs(Vmin - Vmax) + abs(V0 - Vmin);
    tmax = cycles * deltaV / scan_rate;

    disp('doCV.m: performing cyclic voltammogram');
    sol_CV = VappFunction(sol, 'tri', [V0, Vmax, Vmin, cycles, tmax / cycles], tmax, tpoints, 0);
    disp('doCV.m: complete');

end
