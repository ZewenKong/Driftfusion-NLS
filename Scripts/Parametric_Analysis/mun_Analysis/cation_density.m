%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files/pmpi_v2.csv';
cat_density = [1e18, 5e18, 1e19, 5e19, 1e20];
var = cat_density;

par = pc(input);
par.prob_distro_function = 'Boltz';

% - - - - - - - - - - handles
soleqs = cell(1, length(var));
sols = cell(size(soleqs));
maxFs = cell(size(sols));
muns = cell(size(sols));

%% - - - - - - - - - - DATA PROCESSING - - - - - - - - - -

for i = 1:length(var)
    par.Ncat(:) = var(i);
    par = refresh_device(par);
    soleqs{i} = equilibrate(par);
end

%% - - - - - - - - - - DO MEASUREMENTS - - - - - - - - - -

for i = 1:length(soleqs)
    soleq = soleqs{i};
    sol = doCV(soleq.ion, 0, 0, -1, 1, 1e-1, 1, 500); % cycle = 1
    sols{i} = sol;
end

%% - - - - - - - - - - DO CALCULATION - - - - - - - - - -

for i = 1:length(sols)
    sol = sols{i};
    xmesh = sol.x;
    xpos = 0;
    ppos = getpointpos(xpos, xmesh);

    [F_idx, maxF] = dfana_ionic.maxFatIF(sol);

    maxFs{i} = maxF; % save the max field value
end

for i = 1:length(maxFs) % Loop
    maxF = maxFs{i};
    mun = dfcalc_transportmodel.normdmu(maxF);
    muns{i} = mun; % save the normalised mobility
end

%% - - - - - - - - - - PLOTTING - - - - - - - - - -

figure('Name', 'Normalised Mobility v.s. Cation Density');

% - - - - - - - - - - scatter plot
mun_array = cell2mat(muns);
semilogx(cat_density, mun_array, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;
semilogx(cat_density, mun_array, '--', 'LineWidth', 1);
hold off;

xlabel('Cation Density (Ncat)');
ylabel('Normalized Mobility \mu_n');
title('\mu_n v.s. Ncat');
