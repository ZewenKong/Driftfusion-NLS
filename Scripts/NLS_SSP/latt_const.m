%% - - - - - - - - - - CODE START - - - - - - - - - -

input = 'Input_files/pmpi_v2.csv';
lattice = [4.32e-8, 5.32e-8, 6.32e-8, 7.32e-8, 8.32e-8];
var = lattice;

par = pc(input);
par.prob_distro_function = 'Boltz';

soleq = equilibrate(par);
soleq.el.par.isEquilibrate = "sim";
soleq.ion.par.isEquilibrate = "sim";

sol = doCV(soleq.ion, 0, 0, -1, 1, 5e-1, 1, 500);

xmesh = sol.x;
xpos = 0;
ppos = getpointpos(xpos, xmesh);

[F_idx, maxF] = df_analysis.maxFatIF(sol);

muns = cell(1, length(var));

for i = 1:length(var)
    kB = 1.38e-23; T = 298; d = var(i);
    z = 1; e = 1.6e-19;
    F = abs(maxF);
    Ea = 0.45 * e; v0 = 5e12; c = 1;

    P = c * v0 * exp(- (Ea / (kB * T)));

    vD = d * P * (exp((z * e * d * F) / (2 * kB * T)) ...
        - exp(- (z * e * d * F) / (2 * kB * T)));
    vD0 = d * P * ((z * e * d * F) / (kB * T));

    mu = vD ./ F;
    mu0 = vD0 ./ F;
    mun = mu ./ mu0;
    muns{i} = mun;
end

figure('Name', 'Normalised Mobility v.s. Lattice Const');

mun_array = cell2mat(muns);
semilogx(lattice, mun_array, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on
semilogx(lattice, mun_array, '--', 'LineWidth', 1);
hold off

xlabel('Lattice Const. (d)');
ylabel('Normalized Mobility \mu_n');
% title('\mu_n v.s. d');
