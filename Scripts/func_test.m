%% - - - - - - - - - -
% -
% - func. >> ionic charge v.s. time
% - parad. >> dfplot_ionic.Qt_ionic(sol, x1, x2);

dfplot_ionic.Qt_ionic(sol, 0, sol.x(:, end));

%% - - - - - - - - - -
% -
% - func. >> species current at special position point
% - parad. >> dfplot_ionic.Jt_species(sol, position);

dfplot_ionic.Jt_species(sol, sol.x(:, end));

%% - - - - - - - - - -
% -
% - func. >> species current at special time point
% - parad. >> dfplot_ionic.Jt_species(sol, position);
% -
dfplot.Jx(sol);

% ============================================================
% === Energy level at special time points                    =
% ============================================================

%% - - - - - - - - - -
% -
% - func. >> plot eng lvl diagram for special current fluxes
% -
[t_J_max, t_J_corr] = dfana_ionic.tPntOfSpJs(sol, 0);

%% - - - - - - - - - -
% -
% - func. >> energy level diagram at special time point
% -
dfplot.ELx(sol, [t_J_max, t_J_corr]);

%% - - - - - - - - - -
% -
idx_t_J_max = find(t == t_J_max);
idx_t_J_corr = find(t == t_J_corr);

% ============================================================

% ============================================================
% === Special point as new equilibrate solution              =
% ============================================================

%% - - - - - - - - - -
% -
% - func. >> generate equilibrate solution of special points and do CV
% -
[soleq_1, soleq_2] = dfana_ionic.spPntTosoleq(sol, 0);
sol_1 = doCV(soleq_1, 0, 0, -1, 1, 1e-1, 1, 500);
sol_2 = doCV(soleq_2, 0, 0, -1, 1, 1e-1, 1, 500);

%% - - - - - - - - - -
% -
J_1 = dfana.calcJ(sol_1);
Vapp_1 = dfana.calcVapp(sol_1);

xmesh = sol_1.x;
xpos = 0; % define the device position point
ppos_1 = getpointpos(xpos, xmesh);

figure('Name', 'Special point as equilibrium solution')
plot(Vapp_1, J_1.tot(:, ppos_1), 'LineWidth', 1);

%% - - - - - - - - - -
% -
dfplot.ELx(sol_1); % energy level plot
dfplot_ionic.Qt_ionic(sol_1, 0, sol_1.x(:, end)); % charge

%% - - - - - - - - - -
% -
J_2 = dfana.calcJ(sol_2);
Vapp_2 = dfana.calcVapp(sol_2);

xmesh = sol_2.x;
xpos = 0; % define the device position point
ppos_2 = getpointpos(xpos, xmesh);

figure('Name', 'Special point as equilibrium solution')
plot(Vapp_2, J_2.tot(:, ppos_2), 'LineWidth', 1);

%% - - - - - - - - - -
% -
dfplot.ELx(sol_2);
dfplot_ionic.Qt_ionic(sol_2, 0, sol_2.x(:, end));

% ============================================================

% one layer with both silver electrode (w/ or w/o BV)
% const Vset and V read pulse
