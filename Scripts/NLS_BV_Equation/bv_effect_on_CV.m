%% - - - - - - - - - - CODE START - - - - - - - - - -

% read in data files
input = 'Input_files/mapi.csv';
parBV = pc(input);
parBV = refresh_device(parBV); % refresh the device
% parameters setting
cycle = 1; % cycle value
xpos = 0; % define the device position point
area_coeff = 1;
use_abs = 1;
% do equilibrate (BV)
soleqBV = equilibrate(parBV);
soleqBV.el.par.isEquilibrate = "sim";
soleqBV.ion.par.isEquilibrate = "sim";
% do measurements (BV)
solBV = doCV(soleqBV.ion, 0, 0, -1, 1, 0.1, cycle, 500);

% do equilibrate
par = pc(input);
par.j0 = 0; % butler-volmer
par.mobseti = 0; % ionic mobility
par = refresh_device(par);
soleq = equilibrate(par);
soleq.el.par.isEquilibrate = "sim";
soleq.ion.par.isEquilibrate = "sim";
soleq.el.par.j0 = 0;
soleq.ion.par.j0 = 0;
soleq.el.par.mobseti = 0;
soleq.ion.par.mobseti = 0;

% do measurements
sol = doCV(soleq.ion, 0, 0, -1, 1, 0.2, cycle, 500);

%% plot
sols = cell(1, 2);

sols{1} = solBV;
sols{2} = sol;

figure('Name', 'd2d');
hold on;

for i = 1:length(sols)

    J = df_analysis.calcJ(sols{i});
    Vapp = df_analysis.calcVapp(sols{i});
    t = sols{i}.t;
    xmesh = sols{i}.x;
    ppos = getpointpos(xpos, xmesh);

    plot(Vapp, abs(J.tot(:, ppos)), 'DisplayName', ['j0=' num2str(sols{i}.par.j0)], 'LineWidth', 0.5);
end

hold off;
xlabel('Applied Voltage, Vapp [V]');
ylabel('Current Density, J [A cm^{-2}]');
legend('Location', 'northwest', 'FontSize', 10);
yscale log; % logrithm
