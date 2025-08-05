function solstruct = dfII(varargin)
    %
    % adapted from the df.m file,
    % perform the Nernst and Butler-Volmer calculation,
    % to simulate the ionic flux in Driftfusion
    %
    % - - - - - - - - - - CODE START - - - - - - - - - -

    global calc; % check is the 'calc' existing or not,
    clear calc; % clear the 'calc' data structure in workshop
    df_analysis.process_calcdata('reset'); % reset the 'calc' data structure

    % input process
    if length(varargin) == 0
        par = pc;
        dficAnalytical = true;

    elseif length(varargin) == 1
        icsol = varargin{1, 1}.u;
        icx = varargin{1, 1}.x;
        par = varargin{1, 1}.par;
        dficAnalytical = false;

    elseif length(varargin) == 2

        if max(max(max(varargin{1, 1}.u))) == 0
            par = varargin{2};
            dficAnalytical = true;

        elseif isa(varargin{2}, 'char') == 1
            input_solstruct = varargin{1, 1};
            icsol = input_solstruct.u;
            icx = input_solstruct.x;
            par = input_solstruct.par;
            dficAnalytical = false;

        else
            input_solstruct = varargin{1, 1};
            icsol = input_solstruct.u;
            icx = input_solstruct.x;
            par = varargin{2};
            dficAnalytical = false;

        end

    end

    % - - - - - - - - - - UNPACK PROPERTIES - - - - - - - - - -

    % - - - - - driftfusionNLS: model switch [dfII.m]
    isEquilibrate = par.isEquilibrate; % equilibrate 'eq' or simulation 'sim'
    isECM = par.isECM; % electrochemical switch
    k0_trap = par.k0_trap;
    dynamic_adp = par.dynamic_adp;
    % - - - - - END

    % - - - - - driftfusionNLS: physical consts (nernst & butler-volmer) [dfII.m]
    E_st = par.E_st; % standard potential for (Ag + I- <--> AgI + e-, 0.152 eV)
    E_hyd = par.E_hyd; % standard hydrogen electrode (SHE)
    z = par.z; % no. of electrons involved in the electrode reaction
    j0 = par.j0; % exchange current density (right side)
    sa_r = par.sa_r; % ionic charge recombination rate
    alpha_e = par.alpha_c; % electrode charge transfer coefficient
    % - - - - - END

    % temperature
    T = par.T;

    % physical constants
    kB = par.kB;
    q = par.q;
    e = par.e;
    epp0 = par.epp0;
    R = par.R; % universal gas constant
    F = par.F; % Faraday constant
    A = par.A; % Avogadro constant

    % dependent properties
    Vbi = par.Vbi; % built-in voltage
    n0_l = par.n0_l; n0_r = par.n0_r; % equilibrium electron density
    p0_l = par.p0_l; p0_r = par.p0_r; % equilibrium hole density

    dev = par.dev; % device parameters

    % spatial mesh
    xmesh = par.xx; % spatial mesh, thickness of the device
    x_sub = par.x_sub; % spatial mesh (not start from 0)
    x = xmesh;

    % time mesh
    t = meshgen_t(par);

    % device parameters
    N_ionic_species = par.N_ionic_species; % no. of ionic species in this solution (2)
    N_variables = par.N_ionic_species + 3; % no. of variables in this solution (+3 for V, n, and p)
    N_max_variables = par.N_max_variables; % maximum number of variables in this version

    device = par.dev_sub; % device sub parameters

    % charge parameters
    mu_n = device.mu_n; % electron mobility
    mu_p = device.mu_p; % hole mobility
    mu_c = device.mu_c; % cation mobility
    mu_a = device.mu_a; % anion mobility

    Nc = device.Nc; % conduction band eDOS
    Nv = device.Nv; % valence band eDOS

    c_max = device.c_max; % cation density upper limit
    a_max = device.a_max; % anion density upper limit

    gradNc = device.gradNc; % conduction band effective density of states gradient
    gradNv = device.gradNv; % valence band effective density of states gradient
    gradEA = device.gradEA; % electron affinity gradient
    gradIP = device.gradIP; % ionisation potential gradient

    epp = device.epp; % dielectric constant
    epp_factor = par.epp_factor; % maximum dielectric constant (for normalisation)

    ni = device.ni; % intrinsic carrier density

    taun = device.taun; % electron SRH time constant
    taup = device.taup; % hole SRH time constant
    taun_vsr = device.taun_vsr; % electron SRH time constant
    taup_vsr = device.taup_vsr; % hole SRH time constant

    nt = device.nt; % SRH electron trap constant
    pt = device.pt; % SRH hole trap constant

    NA = device.NA; % acceptor doping density
    ND = device.ND; % donor doping density

    switch N_ionic_species

        case 0 % Nani, Ncat, a, and c set to zero for Poisson
            Ncat = zeros(1, length(x_sub));
            Nani = zeros(1, length(x_sub));

        case 1 % Nani and a both set to zero for Poisson
            Ncat = device.Ncat;
            Nani = zeros(1, length(x_sub));

        case 2
            Ncat = device.Ncat;
            Nani = device.Nani;
    end

    xprime_n = device.xprime_n; % translated x co-ordinates for interfaces
    xprime_p = device.xprime_p;

    sign_xn = device.sign_xn; % 1 if xn increasing, -1 if decreasing wrt x
    sign_xp = device.sign_xp; % 1 if xp increasing, -1 if decreasing wrt x

    alpha0_xn = device.alpha0_xn; % alpha0_xn is alpha for F = 0 reference to xprime_n
    beta0_xp = device.beta0_xp; % beta0_xp is beta for F = 0 referenced to xprime_p

    z_c = par.z_c; % +1
    z_a = par.z_a; % -1

    n0_l = par.n0_l; n0_r = par.n0_r; % equilibrium charge (electron, hole) density
    p0_l = par.p0_l; p0_r = par.p0_r;

    sn_l = par.sn_l; sn_r = par.sn_r; % charge surface recombination rate at interface
    sp_l = par.sp_l; sp_r = par.sp_r;

    Rs = par.Rs; % series resistance
    gamma = par.gamma;

    % - - - - - driftfusionNLS: recombination coefficient [dfII.m]
    B = device.B; % electronic
    B_ionic = device.B_ionic; % ionic
    % - - - - - END

    % switches and accelerator coefficients
    mobset = par.mobset; % electronic carrier transport switch
    mobseti = par.mobseti; % ionic carrier transport switch

    K_c = par.K_c; % cation transport rate multiplier
    K_a = par.K_a; % anion transport rate multiplier

    radset = par.radset; % radiative recombination switch (1)
    SRHset = par.SRHset; % SRH recombination switch (1)

    vsr_zone = device.vsr_zone;
    srh_zone = device.srh_zone;

    Rs_initial = par.Rs_initial;
    Field_switch = dev.Field_switch;

    % - - - - - - - - - - GENERATION FUNCTION - - - - - - - - - -

    g1_fun = fun_gen(par.g1_fun_type); % constant
    g2_fun = fun_gen(par.g2_fun_type);

    gxt1 = 0;
    gxt2 = 0;
    g = 0;
    gx1 = par.gx1; % light source 1 & 2
    gx2 = par.gx2;

    int1 = par.int1;
    int2 = par.int2;

    g1_fun_type = par.g1_fun_type; % constant
    g2_fun_type = par.g2_fun_type;

    g1_fun_arg = par.g1_fun_arg; % 0
    g2_fun_arg = par.g2_fun_arg;

    if strcmp(g1_fun_type, 'constant') % string comparison function
        % illumination type g1_fun_type and g2_fun_type convert to Boolean
        % for faster execution in PDEPE
        g1_fun_type_constant = 1;
    else
        g1_fun_type_constant = 0;
    end

    if strcmp(g2_fun_type, 'constant')
        g2_fun_type_constant = 1;
    else
        g2_fun_type_constant = 0;
    end

    gM = g1_fun(g1_fun_arg, t') * gx1 + g2_fun(g2_fun_arg, t') * gx2; % check for negative generation and deal error if present

    if any(any(gM < 0))
        error('dfionic.m: generation cannot be negative - please check your generation function and associated inputs')
    end

    % voltage function
    Vapp_fun = fun_gen(par.V_fun_type);
    Vres = 0;
    J = 0;

    % solver variables
    i = 1;
    V = 0; n = 0; p = 0; a = 0; c = 0;
    dVdx = 0; dndx = 0; dpdx = 0; dadx = 0; dcdx = 0;
    F_V = 0; F_n = 0; F_p = 0; F_c = 0; F_a = 0;
    S_V = 0; S_n = 0; S_p = 0; S_c = 0; S_a = 0;
    r_rad = 0; r_srh = 0; r_vsr = 0; r_np = 0;
    alpha = 0; beta = 0;
    G_n = 1; G_p = 1; % diffusion enhancement prefactor of electrons/holes

    % initialise solution arrays
    u_maxvar = zeros(N_max_variables, 1); % create a zeros matrix with size (N_max_variables, 1), a column vector
    dudx_maxvar = zeros(N_max_variables, 1);
    ul_maxvar = zeros(N_max_variables, 1);
    ur_maxvar = zeros(N_max_variables, 1);

    % - - - - - - - - - - SOLVER OPTIONS - - - - - - - - - -

    options = odeset('MaxStep', par.MaxStepFactor * 0.1 * par.tmax, ... % MaxStep = limit maximum time step size during integration
        'RelTol', par.RelTol, ... % increase the limit to achieve fast pde convergence
        'AbsTol', par.AbsTol); % increase the limit to achieve fast pde convergence

    % - - - - - - - - - - CALL SOLVER - - - - - - - - - -

    % - - - - - driftfusionNLS: model disp [dfII.m]
    disp("dfII.m - model: " + isEquilibrate);
    disp('-');
    % - - - - - END

    % inputs with '@' are function handles to the subfunctions
    % below for the: equation, initial conditions, boundary conditions
    % u, the solution matrix, a 3D matrix for which the dimensions are [time, spacem variables]
    % u icludes, V, n, p, c, a (in order)

    u = pdepe(par.m, @dfpde, @dfic, @dfbc, x, t, options);

    % - - - - - - - - - - OUTPUTS - - - - - - - - - -

    solstruct.u = u; % save 'u' to the 'solstruct.u' structural variable
    solstruct.x = x;
    solstruct.t = t;
    solstruct.par = par; % store parameters object

    if par.vsr_mode == 1 && par.vsr_check == 1
        compare_rec_flux(solstruct, par.RelTol_vsr, par.AbsTol_vsr, 0);
    end

    % - - - - - - - - - - SUBFUNCTIONS - - - - - - - - - -

    % set up partial differential equation (pdepe) (see MATLAB pdepe help for details of C,F,S),
    % C = Time-dependence prefactor; F = Flux terms; S = Source terms; dudx is the MATLAB-created variable

    % - - - - - driftfusionNLS: model disp [dfII.m]
    disp("dfII.m - ECM switch state: " + isECM);
    disp('-');
    % - - - - - END

    function [C, F, S] = dfpde(x, t, u, dudx)

        if x == x_sub(1) % reset position point
            i = 1;
        end

        if g1_fun_type_constant % generation function (illumination)
            gxt1 = int1 * gx1(i);
        else
            gxt1 = g1_fun(g1_fun_arg, t) * gx1(i);
        end

        if g2_fun_type_constant
            gxt2 = int2 * gx2(i);
        else
            gxt2 = g2_fun(g2_fun_arg, t) * gx2(i);
        end

        g = gxt1 + gxt2;

        % unpack variables, assign first 'N_variables' value of u and dudx
        u_maxvar(1:N_variables) = u;
        dudx_maxvar(1:N_variables) = dudx;

        V = u_maxvar(1);
        n = u_maxvar(2);
        p = u_maxvar(3);
        c = u_maxvar(4);
        a = u_maxvar(5);

        dVdx = dudx_maxvar(1); % dVdx, 电场的局部变化率
        dndx = dudx_maxvar(2);
        dpdx = dudx_maxvar(3);
        dcdx = dudx_maxvar(4);
        dadx = dudx_maxvar(5);

        G_n = Nc(i) / (Nc(i) - gamma * n); % diffusion enhancement prefactors (gamma = 0 for Boltz)
        G_p = Nv(i) / (Nv(i) - gamma * p);

        C_V = 0;
        C_n = 1;
        C_p = 1;
        C_c = 1;
        C_a = 1;
        C = [C_V; C_n; C_p; C_c; C_a];

        F_V = (epp(i) / epp_factor) * dVdx;
        % electronic flux term
        F_n = mu_n(i) * n * (-dVdx + gradEA(i)) + (G_n * mu_n(i) * kB * T * (dndx - ((n / Nc(i)) * gradNc(i))));
        F_p = mu_p(i) * p * (dVdx - gradIP(i)) + (G_p * mu_p(i) * kB * T * (dpdx - ((p / Nv(i)) * gradNv(i))));
        % ionic flux term
        F_c = mu_c(i) * (z_c * c * dVdx + kB * T * (dcdx + (c * (dcdx / (c_max(i) - c)))));
        F_a = mu_a(i) * (z_a * a * dVdx + kB * T * (dadx + (a * (dadx / (a_max(i) - a)))));

        F = [F_V; mobset * F_n; mobset * F_p; mobseti * K_c * F_c; mobseti * K_a * F_a];

        r_rad = radset * B(i) * (n * p - ni(i) ^ 2); % radiative
        r_srh = SRHset * srh_zone(i) * ((n * p - ni(i) ^ 2) / (taun(i) * (p + pt(i)) + taup(i) * (n + nt(i)))); % bulk SRH

        alpha = (sign_xn(i) * q * dVdx / (kB * T)) + alpha0_xn(i);
        beta = (sign_xp(i) * q * -dVdx / (kB * T)) + beta0_xp(i);

        r_vsr = SRHset * vsr_zone(i) * ((n * exp(-alpha * xprime_n(i)) * p * exp(-beta * xprime_p(i)) - nt(i) * pt(i)) ...
            / (taun_vsr(i) * (p * exp(-beta * xprime_p(i)) + pt(i)) + taup_vsr(i) * (n * exp(-alpha * xprime_n(i)) + nt(i)))); % volumetric surface recombination

        r_np = r_rad + r_srh + r_vsr; % total electron and hole recombination

        S_V = (1 / (epp_factor * epp0)) * (-n + p - NA(i) + ND(i) + z_a * a + z_c * c - (z_a * Nani(i) + z_c * Ncat(i)));
        S_n = g - r_np;
        S_p = g - r_np;

        % - - - - - driftfusionNLS: ions annihilation and recombination [dfII.m]
        % non-negativity constraint on ion denstiy
        a = max(a, 0); c = max(c, 0);

        % Poole-Frenkel recombination
        delta_ac = a * c - ((dev.Nani(i)) * (dev.Ncat(i)));
        r_iv = radset * B_ionic(i) * (delta_ac);

        if isECM == "ecm_on"
            % ions trapping/detrapping
            ions_bias = abs(a - dev.Nani(i)) / dev.Nani(i) + abs(c - dev.Ncat(i)) / dev.Ncat(i); % anion bias ratio + cation bias ratio
            k_trap = 1 * (k0_trap + (dynamic_adp * ions_bias));
            S_a = (-r_iv) + k_trap * (dev.Nani(i) - a); % (-r_iv), recombination
            S_c = (-r_iv) + k_trap * (dev.Ncat(i) - c); % (k_trap * (dev.Nani(i) - a)), trapping/detrapping
        else
            S_a = (-r_iv);
            S_c = (-r_iv);
        end

        S = [S_V; S_n; S_p; S_c; S_a];
        % - - - - - END

        C = C(1:N_variables); % remove unused variables
        F = F(1:N_variables); %
        S = S(1:N_variables); %

        i = i + 1;
    end

    function u0 = dfic(x)

        % dfic, driftfusion initial condition return the initial condition, u0
        % with data type - a numeric array []

        if x == x_sub(1)
            i = 1;
        end

        if length(par.dcell) == 1 % single layer
            % dcell, numeric array includes cumulative thickness
            u0_ana = [
                      (x / xmesh(end)) * Vbi;
                      n0_l * exp((x * (log(n0_r) - log(n0_l))) / par.dcum0(end)); % electron density, n0(x), changes with the spatial x
                      p0_l * exp((x * (log(p0_r) - log(p0_l))) / par.dcum0(end)); % hole density, p0(x)
                      dev.Ncat(i); % cation density
                      dev.Nani(i); % anion density
                      ];
        else % multi-layered
            u0_ana = [
                      (x / xmesh(end)) * Vbi;
                      dev.n0(i);
                      dev.p0(i);
                      dev.Ncat(i);
                      dev.Nani(i);
                      ];
        end

        u0_ana = u0_ana(1:N_variables); % confirm the u0_ana with size of N_variables

        % organise ICs based on number of variables and SOL_IC
        if dficAnalytical
            u0 = u0_ana;
        else
            u0_input = interp1(icx, squeeze(icsol(end, :, :)), x)'; % initial conditions taken from input solution

            if N_variables > length(u0_input) % if the number of variables has increased then add analytical, ICs for missing variables

                u0(1:length(u0_input), 1) = u0_input; % add initial conditions for new variables from U_ANA
                u0(length(u0_input) + 1:N_variables, 1) = u0_ana(length(u0_input) + 1:N_variables);
            else
                u0 = u0_input;
            end

        end

        i = i + 1;
    end

    function [Pl, Ql, Pr, Qr] = dfbc(xl, ul, xr, ur, t) % Driftfusion boundary condition

        % refer to PDEPE help for the precise meaning of P and Q;
        % l and r refer to left and right boundaries

        ul_maxvar(1:N_variables) = ul;
        ur_maxvar(1:N_variables) = ur;

        V_l = ul_maxvar(1);
        V_r = ur_maxvar(1); % the potential at the right boundary
        n_l = ul_maxvar(2);
        n_r = ur_maxvar(2);
        p_l = ul_maxvar(3);
        p_r = ur_maxvar(3);
        c_l = ul_maxvar(4);
        c_r = ur_maxvar(4); % the cation density at the right boundary
        a_l = ul_maxvar(5);
        a_r = ur_maxvar(5); % the anion density at the right boundary

        % voltage apply function
        switch par.V_fun_type
            case 'constant'
                Vapp = par.V_fun_arg(1);
            otherwise
                Vapp = Vapp_fun(par.V_fun_arg, t); % V_fun_arg = 0
        end

        % (1) Initialisation

        % initial boundary potential (no bias applied)
        % left/right electrode alignment (to SHE, same criterion as E_st)
        E_left_boundary_init = E_hyd - par.Phi_left; % E_hyd = -4.44 (in vaccum), E_st = -0.152 (in SHE)
        E_right_boundary_init = E_left_boundary_init - Vbi; % V_r == Vbi - Vapp - Vres; initially, Vapp is 0, and V_r = Vbi

        % initial equilibrium ion density (Ag+ here)
        a0_r = df_analysis.eqdens(E_right_boundary_init, E_st, R, T, F, A, z);

        % (2) Butler-Volmer electrochemical current/flux calculation

        % (2.1) Nernst equation (the SHE potential)

        % right boundary
        % transfer the density to chemical activity at right boundary
        rho_r = df_analysis.denstochemact(a_r, A);

        % the chemical potential of iodine interstitial at right boundary,
        % used in butler-volmer calculation
        E_right_boundary_eq = df_analysis.nernst(rho_r, E_st, R, T, F, z);

        % left boundary
        rho_l = df_analysis.denstochemact(a_l, A);
        E_left_boundary_eq = df_analysis.nernst(rho_l, E_st, R, T, F, z);

        % (2.2) Butler-Volmer equation

        % in resistive switching, the reaction between iodine interstitial and silver electrode happens,
        % Ag + I- <--> AgI + e-
        % the butler-volmer current arises from the reaction

        % in driftfusion, the potential of the left electrode is set to zero;
        % therefore, an adjustment is added to the standard potential.
        % In details, the left electrode is set to be zero (reference point),
        % but in vaccum the hydrogen is the reference point,
        % the standard potential of the Ag/AgI should be adjusted

        % boundary potential (bias applied)
        E_left_boundary = E_hyd - par.Phi_left; % E_hyd = -4.44 (in vaccum), E_st = -0.152 (in SHE)
        E_right_boundary = E_left_boundary - V_r; % V_r = = Vbi - Vapp - Vres; initially, Vapp is 0, and V_r = Vbi

        % overpotential at right boundary
        eta_r = E_right_boundary - E_right_boundary_eq;

        % butler-volmer current & ionic flux
        j_bv_r = df_analysis.butlervolmer(j0, alpha_e, R, T, F, E_right_boundary, E_right_boundary_eq);

        % save the calculated data
        df_analysis.process_calcdata("t", t, ... % time
            "E_left_boundary_init", E_left_boundary_init, ... % initial left boundary potential
            "E_right_boundary_init", E_right_boundary_init, ... % initial right boundary potential
            "a0_r", a0_r, ... % initial anion density at right boundary
            "a_r", a_r, ... % anion density at right boundary
            "rho_r", rho_r, ... % anion concentration (chemical activity) at right boundary
            "E_right_boundary_eq", E_right_boundary_eq, ... % right boundary potential at equilibrium
            "E_right_boundary", E_right_boundary, ... % right boundary potential
            "eta_r", eta_r, ... % overpotential at right boundary
            "j_bv_r", j_bv_r); % butler-volmer flux at right boundary

        if Rs == 0
            Vres = 0;
        else
            J = e * sp_r * (p_r - p0_r) - e * sn_r * (n_r - n0_r) - j_bv_r; % electron current + hole current + electrochemical (butler-volmer) current

            if Rs_initia
                Vres = -J * Rs * t / par.tmax; % initial linear sweep
            else
                Vres = -J * Rs;
            end

        end

        % check the "dfII.m" is used for equilibrate or simulate
        if isEquilibrate == "eq" % equilibrate
            Pr_a = mobseti * (sa_r * (a_r - a0_r));
        else % simulate (butler-volmer ionic flux)
            Pr_a = (- j_bv_r);
        end

        Pl = [-V_l;
              mobset * (-sn_l * (n_l - n0_l));
              mobset * (-sp_l * (p_l - p0_l));
              0;
              0; ];

        Ql = [0;
              1;
              1;
              1;
              1; ];

        Pr = [-V_r + Vbi - Vapp - Vres; % V_r = Vbi - Vapp - Vres
              mobset * (sn_r * (n_r - n0_r));
              mobset * (sp_r * (p_r - p0_r));
              0;
              Pr_a; ];

        Qr = [0;
              1;
              1;
              1;
              1; ];

        Pl = Pl(1:N_variables); % remove unused entries
        Pr = Pr(1:N_variables);
        Ql = Ql(1:N_variables);
        Qr = Qr(1:N_variables);
    end

end
