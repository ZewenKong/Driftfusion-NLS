%% - - - - - -
% - - - func. >> ionic charge v.s. time
% - - - parad. >> dfplot_ionic.Qt_ionic(sol, x1, x2);
% - - -
dfplot_ionic.Qt_ionic(sol, 0, sol.x(:, end));

%% - - - - - -
% - - - func. >> energy lvl at special time point
% - - - parad. >> dfplot.ELx(sol, time array);
% - - -
dfplot.ELx(sol, [1, 2, 3]);

%% - - - - - -
% - - - func. >> species current at special position point
% - - - parad. >> dfplot_ionic.Jt_species(sol, position);
% - - -
dfplot_ionic.Jt_species(sol, sol.x(:, end));

%% - - - - - -
% - - - func. >> species current at special time point
% - - - parad. >> dfplot_ionic.Jt_species(sol, position);
% - - -
dfplot.Jx(sol);
