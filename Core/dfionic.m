function solstruct = dfionic(varargin)
    %
    % Adapted from the df.m file,
    % perform the Nernst and Butler-Volmer calculation,
    % to simulate the ionic flux.
    %
    %% - - - - - - - - - - CODE START - - - - - - - - - -

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

    %% - - - - - - - - - - UNPACK PROPERTIES - - - - - - - - - -

    % - - - - - - - - - - physics constants
    kB = par.kB; q = par.q; e = par.e; epp0 = par.epp0; T = par.T;

    % - - - - - - - - - - * constants for nernst & butler-volmer calculation
    E_st = par.E_st; R = par.R; F = par.F;
    z = par.z; j0 = par.j0; alpha_e = par.alpha_c;
    E_hyd = par.E_hyd;

    % - - - - - - - - - - spatial mesh
    xmesh = par.xx; % spatial mesh, thickness of the device
    x_sub = par.x_sub; % spatial mesh (not start from 0)
    x = xmesh;

    % - - - - - - - - - - time mesh
    % - - - - - - - - - - original ver.
    t = meshgen_t(par);
    %
    % - - - - - - - - - - update ver. (increase time points (for fast pde convergence))
    % t = meshgen_t(par);
    % par.tpoints = max(1000, par.tpoints * 5);

    % - - - - - - - - - - dependent properties
    Vbi = par.Vbi; % built-in voltage
    n0_l = par.n0_l; n0_r = par.n0_r; % equilibrium electron density
    p0_l = par.p0_l; p0_r = par.p0_r; % equilibrium hole density

    % - - - - - - - - - - device parameters
    N_ionic_species = par.N_ionic_species; % number of ionic species in this solution (2)
    N_variables = par.N_ionic_species + 3; % number of variables in this solution (+3 for V, n, and p)
    N_max_variables = par.N_max_variables; % maximum number of variables in this version

    dev = par.dev; % device parameters
    device = par.dev_sub; % device sub parameters

    % - - - - - - - - - - charge parameters
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

    z_c = par.z_c;
    z_a = par.z_a;

    n0_l = par.n0_l; n0_r = par.n0_r; % equilibrium charge (electron, hole) density
    p0_l = par.p0_l; p0_r = par.p0_r;

    sn_l = par.sn_l; sn_r = par.sn_r; % charge surface recombination rate at interface
    sp_l = par.sp_l; sp_r = par.sp_r;

    Rs = par.Rs; % series resistance
    gamma = par.gamma;

    % - - - - - - - - - - recombination (electronic & ionic)
    B = device.B; % electronic
    Bionic = device.Bionic; % ionic (build_device.m)

    % - - - - - - - - - - switches and accelerator coefficients
    mobset = par.mobset; % electronic carrier transport switch
    mobseti = par.mobseti; % ionic carrier transport switch

    K_c = par.K_c; K_a = par.K_a; % cation/anion transport rate multiplier

    radset = par.radset; % radiative recombination switch (1)
    SRHset = par.SRHset; % SRH recombination switch (1)

    vsr_zone = device.vsr_zone;
    srh_zone = device.srh_zone;

    Rs_initial = par.Rs_initial;
    Field_switch = dev.Field_switch;

    %% - - - - - - - - - - GENERATION FUNCTION - - - - - - - - - -

    g1_fun = fun_gen(par.g1_fun_type); % constant
    g2_fun = fun_gen(par.g2_fun_type);

    gxt1 = 0; gxt2 = 0;
    g = 0;
    gx1 = par.gx1; gx2 = par.gx2; % light source 1 & 2

    int1 = par.int1; int2 = par.int2;

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

    % - - - - - - - - - - voltage function
    Vapp_fun = fun_gen(par.V_fun_type);
    Vres = 0;
    J = 0;

    % - - - - - - - - - - solver variables
    i = 1;
    V = 0; n = 0; p = 0; a = 0; c = 0;
    dVdx = 0; dndx = 0; dpdx = 0; dadx = 0; dcdx = 0;
    F_V = 0; F_n = 0; F_p = 0; F_c = 0; F_a = 0;
    S_V = 0; S_n = 0; S_p = 0; S_c = 0; S_a = 0;
    r_rad = 0; r_srh = 0; r_vsr = 0; r_np = 0;
    alpha = 0; beta = 0;
    G_n = 1; G_p = 1; % diffusion enhancement prefactor of electrons/holes

    % - - - - - - - - - - initialise solution arrays
    u_maxvar = zeros(N_max_variables, 1); % create a zeros matrix with size (N_max_variables, 1), a column vector
    dudx_maxvar = zeros(N_max_variables, 1);
    ul_maxvar = zeros(N_max_variables, 1);
    ur_maxvar = zeros(N_max_variables, 1);

    %% - - - - - - - - - - SOLVER OPTIONS - - - - - - - - - -

    % - - - - - - - - - - original ver.
    % options = odeset('MaxStep', par.MaxStepFactor*0.1*par.tmax, ... % MaxStep = limit maximum time step size during integration
    %     'RelTol', par.RelTol, ...
    %     'AbsTol', par.AbsTol);
    %
    % - - - - - - - - - - * updated ver. (increase the tolerance limit)
    options = odeset('MaxStep', par.MaxStepFactor * 0.1 * par.tmax, ...
        'RelTol', 1e-4, ... % increase the limit to achieve fast pde convergence
        'AbsTol', 1e-7); % increase the limit to achieve fast pde convergence

    %% - - - - - - - - - - CALL SOLVER - - - - - - - - - -

    u = pdepe(par.m, @dfpde, @dfic, @dfbc, x, t, options);

    %% - - - - - - - - - - OUTPUTS - - - - - - - - - -

    solstruct.u = u; % save 'u' to the 'solstruct.u' structural variable
    solstruct.x = x;
    solstruct.t = t;
    solstruct.par = par;

    if par.vsr_mode == 1 && par.vsr_check == 1
        compare_rec_flux(solstruct, par.RelTol_vsr, par.AbsTol_vsr, 0);
    end

    %% - - - - - - - - - - SUBFUNCTIONS - - - - - - - - - -

    function [C, F, S] = dfpde(x, t, u, dudx)

        if x == x_sub(1)
            i = 1;
        end

        if g1_fun_type_constant
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

        u_maxvar(1:N_variables) = u;
        dudx_maxvar(1:N_variables) = dudx;

        V = u_maxvar(1); % 1st variable is V
        n = u_maxvar(2);
        p = u_maxvar(3);
        c = u_maxvar(4);
        a = u_maxvar(5);

        dVdx = dudx_maxvar(1);
        dndx = dudx_maxvar(2);
        dpdx = dudx_maxvar(3);
        dcdx = dudx_maxvar(4);
        dadx = dudx_maxvar(5);

        G_n = Nc(i) / (Nc(i) - gamma * n);
        G_p = Nv(i) / (Nv(i) - gamma * p);

        C_V = 0;
        C_n = 1;
        C_p = 1;
        C_c = 1;
        C_a = 1;
        C = [C_V; C_n; C_p; C_c; C_a];

        F_V = (epp(i) / epp_factor) * dVdx;
        F_n = mu_n(i) * n * (-dVdx + gradEA(i)) + (G_n * mu_n(i) * kB * T * (dndx - ((n / Nc(i)) * gradNc(i))));
        F_p = mu_p(i) * p * (dVdx - gradIP(i)) + (G_p * mu_p(i) * kB * T * (dpdx - ((p / Nv(i)) * gradNv(i))));
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

        r_iv = radset * Bionic(i) * (a * c - ((dev.Nani(i)) * (dev.Ncat(i)))); % radiative

        S_c =- r_iv;
        S_a =- r_iv;

        S = [S_V; S_n; S_p; S_c; S_a];

        C = C(1:N_variables); % remove unused variables
        F = F(1:N_variables); %
        S = S(1:N_variables); %

        i = i + 1;
    end

    function u0 = dfic(x)

        if x == x_sub(1)
            i = 1;
        end

        if length(par.dcell) == 1 % single layer
            u0_ana = [
                      (x / xmesh(end)) * Vbi;
                      n0_l * exp((x * (log(n0_r) - log(n0_l))) / par.dcum0(end)); % electron density, n0(x), changes with the spatial x
                      p0_l * exp((x * (log(p0_r) - log(p0_l))) / par.dcum0(end)); % hole density, p0(x)
                      dev.Ncat(i); % cation density
                      dev.Nani(i); % anion density
                      ];
        else
            u0_ana = [
                      (x / xmesh(end)) * Vbi;
                      dev.n0(i);
                      dev.p0(i);
                      dev.Ncat(i);
                      dev.Nani(i);
                      ];
        end

        u0_ana = u0_ana(1:N_variables);

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

    function [Pl, Ql, Pr, Qr] = dfbc(xl, ul, xr, ur, t)
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

        % - - - - - - - - - - voltage apply function
        switch par.V_fun_type
            case 'constant'
                Vapp = par.V_fun_arg(1);
            otherwise
                Vapp = Vapp_fun(par.V_fun_arg, t); % V_fun_arg = 0
        end

        % - - - - - - - - - - nernst (the SHE potential is used)
        a_boundary = a_r; % the anion density at boundary
        rho_boundary = dfcalc.denstochemact(a_boundary);
        E_eq = dfcalc.nernst(rho_boundary, E_st, R, T, F, z); % SHE potential

        % - - - - - - - - - - butler-volmer electrochemical current (at right hand interface)
        %
        % - - - - - Note
        % In resistive switching, the reaction between iodine interstitial and silver electrode happens,
        % Ag + I- <--> AgI + e-
        % the butler-volmer current arises from the reaction.
        %
        % In driftfusion, the potential of the left electrode is set to zero;
        % therefore, an adjustment is added to the standard potential.
        % In details, the left electrode is set to be zero (reference point),
        % but in vaccum the hydrogen is the reference point,
        % the standard potential of the Ag/AgI should be adjusted.

        E_st_bv = par.Phi_left - E_hyd + E_st;
        E_boundary = E_st_bv + V_r;
        eta = E_boundary - E_eq;

        j_bv = dfcalc.butlervolmer(j0, alpha_e, R, T, F, E_boundary, E_eq); % butler-volmer current density
        f_bv = j_bv / e; % butler-volmer ionic flux

        % fprintf('DEBUG: t=%g, E_boundary=%g, E_eq=%g, f_bv=%g, V_r=%g a_boundary=%g\n Vapp=%g\n', t, E_boundary, E_eq, f_bv, V_r, a_boundary, Vapp);
        % fprintf('DEBUG: t=%g, E_eq_SHE=%g, E_eq_vac=%g, f_bv=%g, V_r=%g a_r=%g\n Vapp=%g\n', t, E_eq_SHE, E_eq_vac, f_bv, V_r, a_r, Vapp);

        % - - - - - - - - - - call a new function to save the calculated data
        dfcalc.save_calcdata(E_boundary, E_eq, eta, a_boundary, t);

        if Rs == 0
            Vres = 0;
        else
            J = e * sp_r * (p_r - p0_r) - e * sn_r * (n_r - n0_r) + j_bv; % electron current +  hole current + electrochemical current

            if Rs_initia
                Vres = -J * Rs * t / par.tmax; % initial linear sweep
            else
                Vres = -J * Rs;
            end

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
              f_bv; ]; % butler-volmer ionic flux

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
