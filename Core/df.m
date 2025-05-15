
function solstruct = df(varargin)
%
% The core DRIFTFUSION function organises properties and inputs for pdepe,
% a routine to test solving the diffusion and drift equations using the matlab pepde solver.
%
% Solution outputs
% V = u(1) = electrostatic potential
% n = u(2) = electron density
% p = u(3) = holes density
% c = u(4) = cation density (optional)
% a = u(5) = anion density (optional)
%
%% - - - - - - - - - - CODE START - - - - - - - - - -

if length(varargin) == 0 % If no input parameter set then call pc directly
    
    par = pc;
    dficAnalytical = true;

elseif length(varargin) == 1 % if one input argument then assume it is the Initial Conditions (IC) solution
    
    icsol = varargin{1, 1}.u; 
    icx = varargin{1, 1}.x;
    par = varargin{1, 1}.par; 
    dficAnalytical = false;

elseif length(varargin) == 2

    if max(max(max(varargin{1, 1}.u))) == 0 % if sol == 0, 'initial solution w/ zero mobility' in 'equilibrate.m'
        
        par = varargin{2}; 
        dficAnalytical = true;
        
    elseif isa(varargin{2}, 'char') == 1 % checks to see if argument is a character
        
        input_solstruct = varargin{1, 1};
        icsol = input_solstruct.u; 
        icx = input_solstruct.x;
        par = input_solstruct.par; 
        dficAnalytical = false;

    else % 'initial solution w/ mobility' and 'initial solution w/ ion mobility' in 'equilibrate.m'
        
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

% - - - - - - - - - - spatial mesh
% xmesh, the spatial mesh
% xx, thickness of the device (without electrode), distribute in number of layer_points
xmesh = par.xx;
% x_sub, thickness of the device (without electrode, 
% not start from zero), distribute in number of layer_points
x_sub = par.x_sub;
x = xmesh;


% - - - - - - - - - - time mesh
% meshgen_t, core function which is used to generate time mesh (t)
t = meshgen_t(par);

% - - - - - - - - - - dependent properties
Vbi = par.Vbi; % built-in voltage
n0_l = par.n0_l; n0_r = par.n0_r; % equilibrium electron density
p0_l = par.p0_l; p0_r = par.p0_r; % equilibrium hole density

% - - - - - - - - - - device parameters
N_ionic_species = par.N_ionic_species; % number of ionic species in this solution (2)
N_variables = par.N_ionic_species + 3; % number of variables in this solution (+3 for V, n, and p)
N_max_variables = par.N_max_variables; % maximum number of variables in this version

dev = par.dev; % device parameters
device = par.dev_sub; % sub-device parameters

% - - - - - - - - - - charge carriers mobilities
mu_n = device.mu_n; % electron mobility
mu_p = device.mu_p; % hole mobility
mu_c = device.mu_c; % cation mobility (e.g. 1e-8 in active layer, define from .csv file)
mu_a = device.mu_a; % anion mobility (e.g. 1e-7 in active layer)

Nc = device.Nc; % conduction band effective density of states (eDOS)
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
taun_vsr = device.taun_vsr; % electron SRH time constant - volumetric interfacial surface recombination scheme
taup_vsr = device.taup_vsr; % hole SRH time constant - volumetric interfacial surface recombination scheme

nt = device.nt; % SRH electron trap constant
pt = device.pt; % SRH hole trap constant

NA = device.NA; % acceptor doping density
ND = device.ND; % donor doping density

switch N_ionic_species

    case 0 % Nani, Ncat, a, and c set to zero for Poisson
        Ncat = zeros(1, length(x_sub)); Nani = zeros(1, length(x_sub));

    case 1 % Nani and a both set to zero for Poisson
        Ncat = device.Ncat; Nani = zeros(1, length(x_sub));

    case 2
        Ncat = device.Ncat; Nani = device.Nani;
end

xprime_n = device.xprime_n; % translated x co-ordinates for interfaces
xprime_p = device.xprime_p;

sign_xn = device.sign_xn; % 1 if xn increasing, -1 if decreasing wrt x
sign_xp = device.sign_xp; % 1 if xp increasing, -1 if decreasing wrt x

alpha0_xn = device.alpha0_xn; % alpha0_xn is alpha for F = 0 reference to xprime_n
beta0_xp = device.beta0_xp; % beta0_xp is beta for F = 0 referenced to xprime_p

z_c = par.z_c; z_a = par.z_a; % +1, -1

n0_l = par.n0_l; n0_r = par.n0_r; % equilibrium charge (electron, hole) density
p0_l = par.p0_l; p0_r = par.p0_r;

sn_l = par.sn_l; sn_r = par.sn_r; % charge surface recombination rate at interface
sp_l = par.sp_l; sp_r = par.sp_r;

Rs = par.Rs; % series resistance
gamma = par.gamma; % blakemore approximation coefficient, 0 for boltzmann stats

% - - - - - - - - - - original ver.
% B = device.B;
%
% - - - - - - - - - - * updated ver.
B = device.B; % radiative recombination rate coefficient
Bionic = device.Bionic; % * rate coefficient for ions and vacancies

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

g1_fun = fun_gen(par.g1_fun_type); % g1_fun_type, is used to control the light source time-dependence (const)
g2_fun = fun_gen(par.g2_fun_type);

gxt1 = 0; gxt2 = 0; g = 0;
gx1 = par.gx1; gx2 = par.gx2; % light source 1 & 2

int1 = par.int1; int2 = par.int2;

g1_fun_type = par.g1_fun_type; % constant
g2_fun_type = par.g2_fun_type;

g1_fun_arg = par.g1_fun_arg; % 0
g2_fun_arg = par.g2_fun_arg;

if strcmp(g1_fun_type, 'constant')

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
    error('df.m: generation cannot be negative - please check your generation function and associated inputs')
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

% - - - - - - - - - - latest ver.
% options = odeset( ...
%     MaxStep=par.MaxStepFactor*0.1*par.tmax, ... % limit maximum time step size during integration
%     RelTol=par.RelTol, ... % relative tolerance (default: 1e-3), controls acceptable relative error
%     AbsTol=par.AbsTol ... % absolute tolerance (default: 1e-6), controls acceptable absolute error
%     );
%
% - - - - - - - - - - * pre-R2021a ver.
options = odeset('MaxStep', par.MaxStepFactor*0.1*par.tmax, ... % MaxStep = limit maximum time step size during integration
    'RelTol', par.RelTol, ...
    'AbsTol', par.AbsTol);

%% - - - - - - - - - - CALL SOLVER - - - - - - - - - -

% inputs with '@' are function handles to the subfunctions
% below for the: equation, initial conditions, boundary conditions
%
% u, the solution matrix, 
% a 3D matrix for which the dimensions are [time, spacem variables]
% u icludes, V, n, p, c, a (in order)

u = pdepe(par.m, @dfpde, @dfic, @dfbc, x, t, options);

%% - - - - - - - - - - OUTPUTS - - - - - - - - - -
% solutions and meshes to structure

solstruct.u = u; % save 'u' to the 'solstruct.u' structural variable
solstruct.x = x;
solstruct.t = t;
solstruct.par = par; % store parameters object

if par.vsr_mode == 1 && par.vsr_check == 1 % volumetric surface recombination error check
    compare_rec_flux(solstruct, par.RelTol_vsr, par.AbsTol_vsr, 0);
end

% - - - - - - - - - - L.J.F.H Code
% if par.vsr_mode == 1 && par.vsr_check == 1
%     try 
%         compare_rec_flux(solstruct, par.RelTol_vsr, par.AbsTol_vsr, 0);
%     catch
%         % Put this here so that probgram doesn't stop for partial solutions
%         % but will stop if soution has failed completely (i.e., only
%         % available for t = 0)
%         if length(solstruct.u(:,1,1)) ~= 1
%             warning("Could not estimate recombination flux error as solution is incomplete")
%         else
%             error("Solution failed, only availabe at t = 0")
%         end
%     end
% end
% - - - - - - - - - -

%% - - - - - - - - - - SUBFUNCTIONS - - - - - - - - - -

% set up partial differential equation (pdepe) (see MATLAB pdepe help for details of C,F,S),
% C = Time-dependence prefactor; F = Flux terms; S = Source terms;
% dudx is the MATLAB-created variable.

    function [C, F, S] = dfpde(x, t, u, dudx)

        % - - - - - - - - - - reset position point
        if x == x_sub(1) % x_sub, the device thickness array
            i = 1;
        end

        % - - - - - - - - - - generation function (illumination)
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

        % - - - - - - - - - - unpack variables
        u_maxvar(1 : N_variables) = u; % N_variables = ionic species (ion, vacancy) + 3 (V, n, p)
                                       % assign first 'N_variables' value of u and dudx to the u_maxvar and dudx_maxvar
        dudx_maxvar(1 : N_variables) = dudx;

        V = u_maxvar(1); % 1st variable is V
        n = u_maxvar(2);
        p = u_maxvar(3);
        c = u_maxvar(4);
        a = u_maxvar(5);

        dVdx = dudx_maxvar(1); % dVdx, 电场的局部变化率 (dVdx, dndx, dpdx, dadx, dcdx are V, n, p, c, a 随 x 的导数)
        dndx = dudx_maxvar(2); % dndx, dpdx, dadx, dcdx, charge carreier density 的局部变化率
        dpdx = dudx_maxvar(3);
        dcdx = dudx_maxvar(4);
        dadx = dudx_maxvar(5);
        
        G_n = Nc(i) / (Nc(i) - gamma * n); % diffusion enhancement prefactors (gamma = 0 for Boltz)
        G_p = Nv(i) / (Nv(i) - gamma * p);
        
        % - - - - - - - - - - equation editor
        % time-dependence pre-factor (pre-allocated above)
        % time-dependence prefactor term
        
        C_V = 0;
        C_n = 1;
        C_p = 1;
        C_c = 1;
        C_a = 1;
        C = [C_V; C_n; C_p; C_c; C_a];
        
        % - - - - - - - - - - flux terms
        %
        % flux, 通量项(物理量在空间上传递或流动的速率)
        
        F_V = (epp(i) / epp_factor) * dVdx;

        % - - - - - Note
        % Electrons flux term
        %
        % First Part: 漂移 Drift
        %
        % mu_n, electron mobility
        % n, electron density
        % (-dVdx + gradEA(i)), 驱动电子运动的驱动力
        % dVdx, 电势梯度产生的电场; gradEA(i), 电子亲和能的梯度用于修正电场
        %
        % Second Part: 扩散 Diffusion
        %
        % G_n, diffusion enhancement prefactor
        % kB * T, thermal energy
        % dndx, 电子密度的空间梯度
        % 在标准的扩散模型中, 电子扩散通量通常与浓度梯度有关, 即可用 dndx 表示.
        % But, 当导带有效态密度 Nc 在空间上不是均匀的时, 电子的浓度 n 与 Nc 之间的归一化关系（n/Nc）会发生变化.
        % ((n/Nc(i)) * gradNc(i))) 的作用是补偿由于 Nc 非均匀变化带来的额外贡献, 梯度修正项.
        % n/Nc(i) 表示在位置 i 处电子浓度相对于该点导带有效态密度的归一化值 (normalisation value).

        F_n = mu_n(i) * n * (-dVdx + gradEA(i)) + (G_n * mu_n(i) * kB * T * (dndx - ((n / Nc(i)) * gradNc(i))));
        F_p = mu_p(i) * p * (dVdx - gradIP(i)) + (G_p * mu_p(i) * kB * T * (dpdx - ((p / Nv(i)) * gradNv(i))));

        % - - - - - Note
        % Cation flux term
        %
        % First Part: 漂移 Drift (z_c * c * dVdx)
        %
        % z_c = 1, cation charge value (+1)
        % c, cation density
        % dVdx, 电势梯度产生的电场
        %
        % Second Part: 扩散 Diffusion
        %
        % kB * T, thermal energy
        % dcdx, 阳离子浓度的梯度 (阳离子浓度 c 关于空间坐标 x 的梯度)
        % (c * (dcdx/(c_max(i) - c))) 项
        % 考虑了当阳离子浓度接近上限 c_max(i) 时, 拥挤效应的修正.
        % 当阳离子浓度 c 增加时, 它们在空间中的排列会变得更加密集; 当浓度接近上限 c_max(i) 时, 离子之间会相互拥挤, 限制它们的自由扩散
        % 乘以 c 可以看作是在高浓度区域中, 拥挤效应对扩散的影响更为显著.
        
        F_c = mu_c(i) * (z_c * c * dVdx + kB * T * (dcdx + (c * (dcdx / (c_max(i) - c)))));
        F_a = mu_a(i) * (z_a * a * dVdx + kB * T * (dadx + (a * (dadx / (a_max(i) - a)))));

        F = [F_V; mobset * F_n; mobset * F_p; mobseti * K_c * F_c; mobseti * K_a * F_a];
        
        % - - - - - - - - - - Electron and hole recombination
        r_rad = radset * B(i) * (n * p - ni(i)^2); % radiative
        r_srh = SRHset * srh_zone(i) * ((n * p - ni(i)^2) / (taun(i) * (p + pt(i)) + taup(i) * (n + nt(i)))); % bulk SRH

        alpha = (sign_xn(i) * q * dVdx / (kB * T)) + alpha0_xn(i);
        beta = (sign_xp(i) * q * -dVdx / (kB * T)) + beta0_xp(i);

        r_vsr = SRHset * vsr_zone(i) * ((n * exp(-alpha * xprime_n(i)) * p * exp(-beta * xprime_p(i)) - nt(i) * pt(i))...
            /(taun_vsr(i) * (p * exp(-beta * xprime_p(i)) + pt(i)) + taup_vsr(i) * (n * exp(-alpha * xprime_n(i)) + nt(i)))); % volumetric surface recombination

        r_np = r_rad + r_srh + r_vsr; % total electron and hole recombination

        % - - - - - - - - - - source terms (V, n, p)
        S_V = (1 / (epp_factor * epp0)) * (-n + p - NA(i) + ND(i) + z_a * a + z_c * c - (z_a * Nani(i) + z_c * Ncat(i)));
        S_n = g - r_np;
        S_p = g - r_np;

        % - - - - - - - - - - Ion and vacancy recombination
        % - - - - - - - - - - original ver. (no ionic recombination)
        % S_c = 0; S_a = 0;
        % 
        % - - - - - - - - - - * updated ver. (ion Frenkel pair recombination)
        % [V_I][I-] <-> ([I_0])^2 % equilibrium reaction at interface PMPbI/PCBM
        r_iv = radset * Bionic(i) * (a * c - ((dev.Nani(i)) * (dev.Ncat(i)))); % radiative 
        S_c = - r_iv;
        S_a = - r_iv;

        S = [S_V; S_n; S_p; S_c; S_a];
        
        C = C(1 : N_variables); % remove unused variables
        F = F(1 : N_variables); %
        S = S(1 : N_variables); %
        
        i = i + 1;
    end

%% - - - - - - - - - - INITIAL CONDITIONS - - - - - - - - - - 

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

        u0_ana = u0_ana(1 : N_variables); % confirm the u0_ana with size of N_variables
        
        % - - - - - - - - - - organise ICs based on number of variables and SOL_IC
        
        if dficAnalytical
            u0 = u0_ana;
        else
            u0_input = interp1(icx, squeeze(icsol(end, :, :)), x)'; % initial conditions taken from input solution
            
            if N_variables > length(u0_input) % if the number of variables has increased then add analytical, ICs for missing variables
                
                u0(1 : length(u0_input), 1) = u0_input; % add initial conditions for new variables from U_ANA
                u0(length(u0_input) + 1 : N_variables, 1) = u0_ana(length(u0_input) + 1 : N_variables);
                
            else
                u0 = u0_input;
            end
        end
        
        i = i + 1;
    end

%% - - - - - - - - - - BOUNDARY CONDITIONS - - - - - - - - - -

% refer to PDEPE help for the precise meaning of P and Q;
% l and r refer to left and right boundaries.

    function [Pl, Ql, Pr, Qr] = dfbc(xl, ul, xr, ur, t)
        
        ul_maxvar(1 : N_variables) = ul;
        ur_maxvar(1 : N_variables) = ur;
        
        V_l = ul_maxvar(1);
        V_r = ur_maxvar(1);
        n_l = ul_maxvar(2);
        n_r = ur_maxvar(2);
        p_l = ul_maxvar(3);
        p_r = ur_maxvar(3);
        c_l = ul_maxvar(4);
        c_r = ur_maxvar(4);
        a_l = ul_maxvar(5);
        a_r = ur_maxvar(5);
        
        switch par.V_fun_type
            case 'constant'
                Vapp = par.V_fun_arg(1);
            otherwise
                Vapp = Vapp_fun(par.V_fun_arg, t);
        end
        
        % flux boundary conditions for both carrier types.
        %
        % calculate series resistance voltage Vres.

        if Rs == 0
            Vres = 0;
        else
            J = e * sp_r * (p_r - p0_r) - e * sn_r * (n_r - n0_r); % here ???
 
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
            0;];
        
        Ql = [0;
            1;
            1;
            1;
            1;];
        
        Pr = [-V_r + Vbi - Vapp - Vres;
            mobset * (sn_r * (n_r - n0_r));
            mobset * (sp_r * (p_r - p0_r));
            0;
            0;];
        
        Qr = [0;
            1;
            1;
            1;
            1;];
        
        Pl = Pl(1 : N_variables); % remove unused entries
        Pr = Pr(1 : N_variables);
        Ql = Ql(1 : N_variables);
        Qr = Qr(1 : N_variables);
    end
end
