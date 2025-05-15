
%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files_zk/peapi.csv';
cat_density = [1e16, 1e17, 1e18, 1e19, 1e20];
var = cat_density;

par = pc(input);
par.prob_distro_function = 'Boltz';
xpos = 0;

% - - - - - - - - - - handles
soleqs = cell(1, length(var));
sols = cell(size(soleqs));
maxFs = cell(size(sols));
muns = cell(size(sols));

%% - - - - - - - - - - DATA PROCESSING - - - - - - - - - -

for i = 1 : length(var)
    par.Ncat(:) = var(i);
    par = refresh_device(par);
    soleqs{i} = equilibrate(par);
end

%% - - - - - - - - - - DO MEASUREMENTS - - - - - - - - - -

for i = 1 : length(soleqs)
    soleq = soleqs{i};
    sol = doCV(soleq.ion, 0, 0, 1, -1, 5e-1, 1, 500); % cycle = 1
    sols{i} = sol;
end

%% - - - - - - - - - - DO CALCULATION - - - - - - - - - -

for i = 1 : length(sols)
    sol = sols{i};
    xmesh = sol.x;
    ppos = getpointpos(xpos, xmesh);
    
    % - - - - - - - - - - Single-layer
    % F = dfana.calcF(sol, "whole");
    % maxF = max(abs(F), [], 'all');
    % - - - - - - - - - -

    % - - - - - - - - - - Three-layer device (interface)
    [F_idx, maxF] = dfcalc_transportmodel.maxFatIF(sol);
    % - - - - - - - - - - 

    % - - - - - - - - - - Three-layer device (within RPP)
    % maxF = maxFwithinRPP(sol); 
    % - - - - - - - - - - 

    maxFs{i} = maxF; % save the max field value
end

for i = 1 : length(maxFs) % Loop
    maxF = maxFs{i};
    mun = dfcalc_transportmodel.normdmu(maxF);
    muns{i}  = mun; % save the normalised mobility
end

%% - - - - - - - - - - Plots

figure('Name', 'Normalised Mobility v.s. Ion Density');

% - - - - - - - - - - scatter plot
mun_array = cell2mat(muns);
semilogx(cat_density, mun_array, 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
hold on;
semilogx(cat_density, mun_array, '--', 'LineWidth', 1);
hold off;

% - - - - - - - - - - fitting curve
% xlog = log10(cat_density);
% y = mun_array;
% p = polyfit(xlog, y, 1);
% xfit = logspace(log10(min(cat_density)), log10(max(cat_density)), 100);
% yfit = polyval(p, log10(xfit));
% 
% hold on
% plot(xfit, yfit, '--', 'LineWidth', 0.5, 'Color', 'black');
% hold off

xlabel('Cation Density (Ncat)');
ylabel('Normalized Mobility \mu_n');
title('\mu_n v.s. Ncat');

current_xlim = xlim; 
x_range = diff(current_xlim); 
increased_xlim = [current_xlim(1) / (1 + 0.5), current_xlim(2) * (1 + 0.5)];
xlim(increased_xlim);

current_ylim = ylim;
y_range = diff(current_ylim);
increased_ylim = [current_ylim(1) / (1 + 0.025), current_ylim(2) * (1 + 0.025)];
ylim(increased_ylim);