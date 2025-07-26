function soleq = equilibrate(varargin)

    % Uses initial conditions defined in DF and runs to equilibrium
    %
    % VARARGIN{1, 1} = PAR
    % VARARGIN{1, 2} = ELECTRONIC_ONLY
    % ELECTRONIC_ONLY:
    % 0 = runs full equilibrate protocol
    % 1 = skips ion equilibration

    %% - - - - - - - - - - CODE START - - - - - - - - - -

    if length(varargin) == 1 % varargin is a cell array

        par = varargin{1, 1}; % returns the first value in varargin (1st row, 1st coloumn) = par
        electronic_only = 0;

    elseif length(varargin) == 2

        par = varargin{1, 1};
        electronic_only = varargin{1, 2}; % returns the second value in varargin (1st row, 2nd coloumn)

    else
        par = pc;
        electronic_only = 0;
    end

    tic; % start stopwatch

    % INITIAL ARGUMENTS & GENERAL PARAMETERS

    % initial arguments

    sol.u = 0;
    % setting sol.u = 0 enables a parameters structure to be read into DF,
    % but indicates that the initial conditions should be the analytical solutions

    par_origin = par; % store the original parameter set

    par.SRHset = 0; % start with zero SRH recombination
    par.radset = 1; % radiative recombination could initially be set to zero in addition if required

    par.N_ionic_species = 0; % start with no ionic carriers
    par.vsr_check = 0; % switch off volumetric surface recombination check

    % general initial parameters

    % voltage apply function & arguments
    par.V_fun_type = 'constant'; % set applied bias to zero
    par.V_fun_arg(1) = 0;

    par.int1 = 0; % set light intensities to zero (source 1 & 2)
    par.int2 = 0;
    par.g1_fun_type = 'constant';
    par.g2_fun_type = 'constant';

    par.tmesh_type = 2; % time mesh
    par.tpoints = 200;

    par.Rs = 0; % series resistance

    % INITIAL SOLUTION W/ ZERO MOBILITY

    % first, the simulation runs with zero mobilities,
    % so that an analytical or initial condition is established
    % without interference from transport effects

    par.mobset = 0; % switch off mobilities (electronic)
    par.mobseti = 0; % switch off mobilities (ionic)

    disp('equilibrate.m: solution initialisation (zero mobility)'); disp('-');
    sol = df(sol, par); % CALL df.m
    disp('equilibrate.m: complete initialisation (zero mobility)'); disp('-');

    % INITIAL SOLUTION W/ MOBILITY

    % by enabling electronic mobilities and activating recombination mechanisms,
    % the system allows electron and hole distributions
    % to adjust and reach an equilibrium state.

    par.mobset = 1; % switch on electronic mobilities
    par.radset = 1; % radiative recombination
    par.SRHset = 1; % shockley–read–hall recombination

    t_diff = (par.dcum0(end) ^ 2) / (2 * par.kB * par.T * min(min(par.mu_n), min(par.mu_p))); % characteristic diffusion time
    par.tmax = 100 * t_diff;
    par.t0 = par.tmax / 1e6;

    disp('equilibrate.m: solution initialisation (electronic mobility)'); disp('-');
    sol = df(sol, par); % CALL df.m

    all_stable = verifyStabilization(sol.u, sol.t, 0.7); % check the solution reached a stabilized status (solution matrix, time array, time increment fraction)

    j = 1;

    while any(all_stable) == 0

        % loop the solution, to check electrons have reached stable config,
        % if not accelerate ions by order of mag,
        % adding up the equilibration time to achieve equilibrium

        disp(['equilibrate.m: increasing equilibration time, tmax = ', num2str(par.tmax * 10 ^ j)]); disp('-');
        par.tmax = 10 * par.tmax;
        par.t0 = par.tmax / 1e6;

        sol = df(sol, par);

        all_stable = verifyStabilization(sol.u, sol.t, 0.7);
    end

    soleq.el = sol; % solution equilibrium (electronic)

    % manually check final section of solution for VSR self-consitency,
    % and compare the interfacial recombination fluxes

    sol_ic = extract_IC(soleq.el, [soleq.el.t(end) * 0.7, soleq.el.t(end)]);
    compare_rec_flux(sol_ic, par.RelTol_vsr, par.AbsTol_vsr, 0);

    soleq.el.par.vsr_check = 1; % switch VSR check on for future use

    disp('equilibrate.m: complete initialisation (electronic mobility)'); disp('-');

    % INITIAL SOLUTION W/ ION MOBILITY

    if electronic_only == 0 && par_origin.N_ionic_species > 0

        % the simulation achieves full electrochemical equilibrium,
        % including both electronic and ionic responses

        par.N_ionic_species = par_origin.N_ionic_species; % par_origin = par
        sol = soleq.el; % create temporary solution for appending initial conditions to
        par.Rs = 0; % par.SRHset = 0; % start without SRH or series resistance

        disp('equilibrate.m: solution initialisation (electronic and ionic mobility)'); disp('-');

        % original ver.
        %
        % only allows for non-zero mobility in active layer),
        % take ratio of electron and ion mobilities in the active layer
        %
        % rat_anion = par.mu_n(par.active_layer)/par.mu_a(par.active_layer);
        % rat_cation = par.mu_n(par.active_layer)/par.mu_c(par.active_layer);
        %
        % updated ver. (accounts non-zero mobility in any layer)
        [max_mu_a, max_mu_a_idx] = max(par.mu_a);
        [max_mu_c, max_mu_c_idx] = max(par.mu_c);
        rat_anion = par.mu_n(par.active_layer) / par.mu_a(max_mu_a_idx);
        rat_cation = par.mu_n(par.active_layer) / par.mu_c(max_mu_c_idx);

        % if the ratio is infinity (ion mobility set to zero),
        % then set the ratio to zero instead

        if isnan(rat_anion) || isinf(rat_anion)
            rat_anion = 0;
        end

        if isnan(rat_cation) || isinf(rat_cation)
            rat_cation = 0;
        end

        par.mobset = 1;
        par.mobseti = 1;
        par.K_a = rat_anion;
        par.K_c = rat_cation;
        par.tmax = 1e4 * t_diff;
        par.t0 = par.tmax / 1e3;

        % original ver.
        % sol = df(sol, par);
        %
        % updated ver.
        sol = dfNLS(sol, par);

        all_stable = verifyStabilization(sol.u, sol.t, 0.7);

        % j0_temp = par.j0; % temporary store of j0 value
        % par.j0 = 0; %set bv j0 to zero for first ion equilibration
        % par = refresh_device(par); % get the value

        % run ionic stabilisation without Butler volmer boundary
        while any(all_stable) == 0

            % loop the solution, to check ions have reached stable config,
            % if not accelerate ions by order of mag

            disp(['equilibrate.m: increasing equilibration time, tmax = ', num2str(par.tmax * 10 ^ j)]); disp('-');

            par.tmax = par.tmax * 10;
            par.t0 = par.tmax / 1e6;

            % original ver.
            % sol = df(sol, par);
            %
            % updated ver.
            sol = dfNLS(sol, par);

            all_stable = verifyStabilization(sol.u, sol.t, 0.7);
        end

        soleq.ion = sol; % write solution
        sol_ic = extract_IC(soleq.ion, [soleq.ion.t(end) * 0.7, soleq.ion.t(end)]); % manually check solution for VSR self-consitency
        compare_rec_flux(sol_ic, par.RelTol_vsr, par.AbsTol_vsr, 0);

        soleq.ion.par.vsr_check = 1; % reset switches
        soleq.ion.par.mobseti = 1;
        soleq.ion.par.K_a = 1;
        soleq.ion.par.K_c = 1;

        disp('equilibrate.m: complete initialisation (electronic and ionic mobility)');
        disp('-');
    end

    disp('equilibrate.m: EQUILIBRATION COMPLETE')
    disp('-');
    toc

end
