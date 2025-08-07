% - - - - - - - - - - CODE START - - - - - - - - - -
% file input
input = "Input_files/single_layer.csv";

% measurement setting
light_intensity = 0;
V0 = 0; V_max = 1; V_min = -1;
scan_rt = 0.2;
cycle = 1;
time_pnts = 500;

%% ECM
par_ecm = pc(input); 
par_ecm.mobseti = 1;
par_ecm.B_ionic = 1e-12;
par_ecm.k0_trap = 1e-6;
par_ecm.dynamic_adp = 1e-9;

par_ecm.RelTol = 1e-6;
par_ecm.AbsTol = 1e-8;
par_ecm = refresh_device(par_ecm);

soleq_ecm = equilibrate(par_ecm);
soleq_ecm.eq.par.isEquilibrate = "sim";
soleq_ecm.ion.par.isEquilibrate = "sim";

sol_ecm = doCV(soleq_ecm.ion, light_intensity, V0, V_max, V_min, scan_rt, cycle, time_pnts);

%% General
par = pc(input);

par.mobseti = 0;
par.B_ionic = 0;
par_ecm.k0_trap = 1e-6;
par_ecm.dynamic_adp = 1e-9;
par.RelTol = 1e-6;
par.AbsTol = 1e-8;
par = refresh_device(par);

soleq = equilibrate(par);
soleq.eq.par.isEquilibrate = "sim";
soleq.ion.par.isEquilibrate = "sim";

sol = doCV(soleq.ion, light_intensity, V0, V_max, V_min, scan_rt, cycle, time_pnts);

%%
% xpos = 0;
% J = df_analysis.calcJ(sol);
% Vapp = df_analysis.calcVapp(sol);
% t = sol.t;
% xmesh = sol.x;
% ppos = getpointpos(xpos, xmesh);
% plot(Vapp, J.tot(:, ppos));
% 


%% Q_ionic v.s. Vapp
% figure('Name', 'Q_ionic v.s. Vapp');
% hold on;
% 
% for i = 1:length(sols)
%     [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sols{i});
%     p1 = find(x <= x1); p1 = p1(end);
%     p2 = find(x <= x2); p2 = p2(end);
%     rho_ionic = df_analysis.rhoIonicCalc(sols{i}, "whole");
%     Vapp = df_analysis.calcVapp(sol);
%     Q_ionic = par.e * trapz(x(p1:p2), rho_ionic(:, p1:p2), 2);
%     plot(Vapp, Q_ionic, 'DisplayName', labels{i});
% end

% xlabel('V_app [V]');
% ylabel('Q_ionic [C cm^{-2}]');
% legend show;
% hold off;

% %% Q_ionic v.s. time
% figure('Name', 'Q_ionic v.s. time');
% hold on;
% 
% for i = 1:length(sols)
%     [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sols{i});
%     rho_ionic = df_analysis.rhoIonicCalc(sols{i}, "whole");
%     p1 = find(x <= x1); p1 = p1(end); % start
%     p2 = find(x <= x2); p2 = p2(end); % end
%     Q_ionic = par.e * trapz(x(p1:p2), rho_ionic(:, p1:p2), 2);
%     plot(t, Q_ionic, 'DisplayName', labels{i});
% end

% xlabel('Time [s]')
% ylabel('Charge [C cm-2]')
% xlim([t(1), t(end)])
% legend('show');
% hold off;

%% J v.s. Vapp

labels = {'sol\_bv', 'sol'};
x1 = 0; x2 = 5e-5; xpos = 0;
sols = cell(1, 2);
sols{1} = sol_ecm;
sols{2} = sol;
figure('Name', 'J v.s. Vapp (BV c.f.)');
hold on;

for i = 1:length(sols)
    J = df_analysis.calcJ(sols{i});
    Vapp = df_analysis.calcVapp(sols{i});
    t = sols{i}.t;
    xmesh = sols{i}.x;
    ppos = getpointpos(xpos, xmesh);
    % plot(Vapp, J.a(:, ppos), 'DisplayName', labels{i});
    % plot(Vapp, J.c(:, ppos), 'DisplayName', labels{i});
    plot(Vapp, J.tot(:, ppos), 'DisplayName', labels{i});
end

hold off;
xlabel('Applied Voltage, Vapp [V]');
ylabel('Current Density, J [A cm^{-2}]');
legend('Location', 'northwest', 'FontSize', 10);

%% ELx
% figure('Name', 'Energy Level Diagram');
% hold on;
% df_plot.ELx_v2(sol_bv);
% df_plot.ELx_v2(sol);
% hold off;

%% npacx
% df_plot.npxacx(sol_bv);
% df_plot.npxacx(sol);

% intrinsic
% no current when bv turns off => resistor
% bv turns on ==> capacitor
