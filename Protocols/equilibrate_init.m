function soleq = equilibrate_init(varargin)
%
% Edit from equilibrate.m
% equilibrate_init() >> initialise
%
%% - - - - - - - - - - CODE START - - - - - - - - - -

if length(varargin) == 1
    par = varargin{1, 1}; electronic_only = 0;
elseif length(varargin) == 2
    par = varargin{1, 1}; electronic_only = varargin{1, 2};
else
    par = pc; electronic_only = 0;
end

tic; % start stopwatch

%% - - - - - - - - - - INITIAL ARGUMENTS & GENERAL PARMEMTERS - - - - - - - - - -

sol.u = 0; par_origin = par;
par.SRHset = 0; par.radset = 1;
par.N_ionic_species = 0;
par.vsr_check = 0;

par.V_fun_type = 'constant'; par.V_fun_arg(1) = 0;
par.int1 = 0; par.int2 = 0; par.g1_fun_type = 'constant'; par.g2_fun_type = 'constant';
par.tmesh_type = 2; par.tpoints = 200;
par.Rs = 0;

%% - - - - - - - - - - INITIAL SOLUTION W/ ZERO MOBILITY - - - - - - - - - -

par.mobset = 0; par.mobseti = 0;

disp('equilibrate_init.m: solution initialisation (zero mobility)'); disp('-');
sol = df(sol, par);
disp('equilibrate_init.m: complete initialisation (zero mobility)'); disp('-');

%% - - - - - - - - - - INITIAL SOLUTION W/ MOBILITY - - - - - - - - - -

par.mobset = 1; par.radset = 1; par.SRHset = 1;
t_diff = (par.dcum0(end)^2) / (2 * par.kB * par.T * min(min(par.mu_n), min(par.mu_p)));
par.tmax = 100 * t_diff; par.t0 = par.tmax / 1e6;

disp('equilibrate_init.m: solution initialisation (electronic mobility)'); disp('-');
sol = df(sol, par);
j = 1;
all_stable = verifyStabilization(sol.u, sol.t, 0.7); 

while any(all_stable) == 0
    disp(['equilibrate_init.m: increasing equilibration time, tmax = ', num2str(par.tmax * 10^j)]); disp('-');
    par.tmax = 10 * par.tmax; par.t0 = par.tmax / 1e6;
    sol = df(sol, par);
    all_stable = verifyStabilization(sol.u, sol.t, 0.7);
end

soleq.el = sol;
sol_ic = extract_IC(soleq.el, [soleq.el.t(end) * 0.7, soleq.el.t(end)]);
compare_rec_flux(sol_ic, par.RelTol_vsr, par.AbsTol_vsr, 0);
soleq.el.par.vsr_check = 1;

disp('equilibrate_init.m: complete initialisation (electronic mobility)'); disp('-');

%% - - - - - - - - - - INITIAL SOLUTION W/ ION MOBILITY - - - - - - - - - -

if electronic_only == 0 && par_origin.N_ionic_species > 0

    par.N_ionic_species = par_origin.N_ionic_species;
    sol = soleq.el;
    par.Rs = 0;

    disp('equilibrate_init.m: solution initialisation (electronic and ionic mobility)'); disp('-');

    % - - - - - - - - - - Original ver. 
    % 
    % Only allows for non-zero mobility in active layer),
    % take ratio of electron and ion mobilities in the active layer.
    % 
    % rat_anion = par.mu_n(par.active_layer)/par.mu_a(par.active_layer);
    % rat_cation = par.mu_n(par.active_layer)/par.mu_c(par.active_layer);
    % - - - - - - - - - - 
    % 
    % - - - - - - - - - - * Updated ver. (accounts non-zero mobility in any layer)
    [max_mu_a, max_mu_a_idx] = max(par.mu_a);
    [max_mu_c, max_mu_c_idx] = max(par.mu_c);

    rat_anion = par.mu_n(par.active_layer) / par.mu_a(max_mu_a_idx);
    rat_cation = par.mu_n(par.active_layer) / par.mu_c(max_mu_c_idx);
    % - - - - - - - - - -

    if isnan(rat_anion) || isinf(rat_anion)
        rat_anion = 0;
    end
    if isnan(rat_cation) || isinf(rat_cation)
        rat_cation = 0;
    end

    par.mobset = 1; par.mobseti = 1;
    par.K_a = rat_anion; par.K_c = rat_cation;
    par.tmax = 1e4 * t_diff; par.t0 = par.tmax / 1e3;

    sol = df(sol, par);

    all_stable = verifyStabilization(sol.u, sol.t, 0.7);

    while any(all_stable) == 0
        disp(['equilibrate_init.m: increasing equilibration time, tmax = ', num2str(par.tmax * 10^j)]); disp('-');
        par.tmax = par.tmax * 10; par.t0 = par.tmax / 1e6;
        sol = df(sol, par);
        all_stable = verifyStabilization(sol.u, sol.t, 0.7);
    end

    soleq.ion = sol;
    sol_ic = extract_IC(soleq.ion, [soleq.ion.t(end) * 0.7, soleq.ion.t(end)]);
    compare_rec_flux(sol_ic, par.RelTol_vsr, par.AbsTol_vsr, 0);
    soleq.ion.par.vsr_check = 1;
    soleq.ion.par.mobseti = 1;
    soleq.ion.par.K_a = 1; soleq.ion.par.K_c = 1;

    disp('equilibrate_init.m: complete initialisation (electronic and ionic mobility)'); disp('-');
end

disp('equilibrate_init.m: EQUILIBRATION COMPLETE'); disp('-');

toc

end