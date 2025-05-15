
%% - - - - - - - - - - CODE START - - - - - - - - - -

% - - - - - - - - - - data inputs
input = 'Input_files_zk/peapi.csv';
cat_density = [1e16, 1e17, 1e18, 1e19, 1e20];
lattice = [5.32e-8, 5.82e-8, 6.32e-8, 6.82e-8, 7.32e-8];
var1 = cat_density;
var2 = lattice;

par = pc(input);
par.prob_distro_function = 'Boltz';
xpos = 0;

%% - - - - - - - - - - VARIABLE 1 - - - - - - - - - -

% - - - - - - - - - - handles
soleqs = cell(1, length(var1));
sols = cell(size(soleqs));
maxFs = cell(size(sols));

% - - - - - - - - - - data processing
for i = 1 : length(var1)
    par.Ncat(:) = var1(i);
    par = refresh_device(par);
    soleqs{i} = equilibrate(par);
end

% - - - - - - - - - - do measurements
for i = 1 : length(soleqs)
    soleq = soleqs{i};
    sol = doCV(soleq.ion, 0, 0, 1, -1, 5e-1, 1, 500); % cycle = 1
    sols{i} = sol;
end

% - - - - - - - - - - do calculation
for i = 1 : length(sols)
    sol = sols{i};
    xmesh = sol.x;
    ppos = getpointpos(xpos, xmesh);

    % - - - - - - - - - - Three-layer device (interface)
    [F_idx, maxF] = dfcalc_transportmodel.maxFatIF(sol);
    % - - - - - - - - - - 

    maxFs{i} = maxF; % save the max field value
end

%% - - - - - - - - - - VARIABLE 2 - - - - - - - - - -

% - - - - - - - - - - handles
muns = cell(length(var1), length(var2));

% - - - - - - - - - - do calculation
for i = 1 : length(var1)

    maxF = maxFs{i};

    for j = 1 : length(var2)

        kB = 1.38e-23; T = 298; z = 1; e = 1.6e-19;
        Ea = 0.45 * e;
        v0 = 5e12; c = 1;

        d = var2(j);
        F = abs(maxF);
        
        P = c * v0 * exp(-(Ea / (kB * T)));

        vD = d * P * (exp((z * e * d * F) / (2 * kB * T)) ...
                - exp(-(z * e * d * F) / (2 * kB * T)));
        vD0 = d * P * ((z * e * d * F)/(kB * T));

        mu = vD ./ F;
        mu0 = vD0 ./ F;
    
        mun = mu ./ mu0;
        muns{i, j} = mun;
    end
end

%% - - - - - - - - - - PLOTTING - - - - - - - - - -

mun_matrix = cell2mat(muns);

log_var1 = log10(var1);
log_var2 = log10(var2);

[X, Y] = meshgrid(log_var1, log_var2);
[Xq, Yq] = meshgrid(...
    linspace(min(log_var1), max(log_var1), 300), ...
    linspace(min(log_var2), max(log_var2), 300));

interp = griddata(X(:), Y(:), mun_matrix(:), Xq, Yq, 'cubic');

figure('Name', 'Focused Contour Plot: Cat & Lat');

imagesc([10^min(log_var1), 10^max(log_var1)], ...
        [10^min(log_var2), 10^max(log_var2)], ...
        interp);

set(gca, 'XScale', 'log', 'YScale', 'log');

colormap(turbo);
colorbar;

clim([1.000, 1.005]);

xlabel('Cation Mobility (cm^2/Vs)');
ylabel('Lattice Const. (cm)');
title('\mu/\mu_0 Contour');
xlim([min(var1), max(var1)]);
ylim([min(var2), max(var2)]);

shading interp;