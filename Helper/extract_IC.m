function sol_ic = extract_IC(sol, requested_time)
    % EXTRACT INITIAL CONDITIONS
    % This function takes a single time point from the solution SOL closest to
    % REQUESTED_TIME and outputs a solution SOL_IC with only that time point. This can
    % then be used as the initial conditions for a new simulation with DRIFTFUSION. The time
    % array SOL_IC.T is set to zero (i.e. a single time point) and VAPP is
    % calculated from the SOL parameters and REQUESTED_TIME.
    % Update 02/09/21 REQUESTED_TIME can now be a 2 element vector [START_TIME,
    % END_TIME]

    if length(requested_time) == 1
        index = 0;
        index = find(sol.t <= requested_time);
        index = index(end);
    elseif length(requested_time) == 2
        index = [0, 0];
        index_temp1 = find(sol.t <= requested_time(1));
        index(1) = index_temp1(end);
        index_temp2 = find(sol.t <= requested_time(2));
        index(2) = index_temp2(end);
    end

    % disp('----- DEBUG extract_IC -----');
    % disp(['size(sol.u)      = [' num2str(size(sol.u)) ']']);
    % disp(['length(sol.t)    = ' num2str(length(sol.t))]);
    % disp(['requested_time   = [' num2str(requested_time) ']']);
    % disp(['computed index   = [' num2str(index) ']']);

    sol_ic = sol;
    sol_ic.u = sol.u(index, :, :); % error

    % Overwrite Vapp
    Vappt = df_analysis.calcVapp(sol);
    sol_ic.par.V_fun_arg(1) = Vappt(index(end));

    % Overwrite time array
    sol_ic.t = sol.t(index);

end
