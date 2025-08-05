%% - - - - - - - - - - CODE START - - - - - - - - - -

input = 'Input_files/pmpi.csv';
% cat_density = [1e16, 1e17, 1e18, 1e19, 1e20];
cat_density = [1e16, 1e18];
var = cat_density;

par = pc(input);
par.prob_distro_function = 'Boltz';

soleqs = cell(1, length(var));
sols = cell(size(soleqs));
maxFs = cell(size(sols));
muns = cell(size(sols));

for i = 1:length(var)
    par.Ncat(:) = var(i);
    par = refresh_device(par);
    soleqs{i} = equilibrate(par);

    soleqs{i}.el.par.isEquilibrate = "sim";
    soleqs{i}.ion.par.isEquilibrate = "sim";
end

for i = 1:length(soleqs)
    soleq = soleqs{i};
    sol = doCV(soleq.ion, 0, 0, -1, 1, 5e-1, 1, 500); % cycle = 1
    sols{i} = sol;
end

for i = 1:length(sols)
    sol = sols{i};
    xmesh = sol.x;
    xpos = 0;
    ppos = getpointpos(xpos, xmesh);

    [F_idx, maxF] = df_analysis.maxFatIF(sol);

    maxFs{i} = maxF; % save the max field value
end

for i = 1:length(maxFs) % Loop
    maxF = maxFs{i};
    mun = df_analysis.normdmu(maxF);
    muns{i} = mun; % save the normalised mobility
end

figure('Name', 'Normalised Mobility v.s. Cation Density');

mun_array = cell2mat(muns);
semilogx(cat_density, mun_array, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;
semilogx(cat_density, mun_array, '--', 'LineWidth', 1);
hold off;

xlabel('Cation Density (Ncat)');
ylabel('Normalized Mobility \mu_n');
title('\mu_n v.s. Ncat');
