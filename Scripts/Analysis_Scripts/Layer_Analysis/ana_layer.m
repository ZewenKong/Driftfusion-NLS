%% - - - - - - - - - - CODE START - - - - - - - - - -

% time point of special current point

% - - - - - - - - - - inputs
xpos = 0;
V_sp = 0.576577; % without BV
% V_sp = 0.28028; % with BV

%% - - - - - - - - - - data processing
t_sp = dfana_ionic.tOfSpV(sol, xpos, V_sp);

%% - - - - - - - - - - energy level at special time point
dfplot.ELxnpxacx(sol, t_sp);
