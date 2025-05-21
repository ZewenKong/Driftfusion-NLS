% The main.m file is responsible for call the function

%% Scan rate confer

%% 1) Compare different scan rate to define the scan rate with the best performance.

scan_rt_cf;

%% 2) I - t curve, to find the best voltage operation range.


%% 3) I - V curve (one cycle), to find the HRS and LRS.

% The I-V curve can also be used to define the ON/OFF ratio and SET/RESET
% voltage.

cyc_to_cyc_scan;

%% 4) Log(I) - Log(V) curve, to find the SET voltage point.


%% 5) I - cycle curve, is used to conduct the endurance test.

% To find the stability of the current during cycles

%% 6) I - V curve (hundreds cycles), is also used to conduct the endurance test



%% 7) I - V curve (D2D)




%% Cycle-to-cycle scan

% I - V Plot
% 
% --> ON/OFF ratio
% --> SET/RESET voltage
% --> Endurance test
% --> Stochasticity/Variability test

%% Device-to-device scan (Different active material)

% I - V Plot
% 
% --> ON/OFF ratio
% --> SET/RESET voltage
% --> Stochasticity/Variability test

dev_to_dev_scan_mat;

%% Device-to-device scan (Different physical properties)

% I - V Plot
% 
% --> ON/OFF ratio
% --> SET/RESET voltage

dev_to_dev_scan_prop;

%% Planned staff

% * Retention time
%
% --> R - t plot, measures long-term retention stability over time
% 
% --> I - t plot, monitors current decay over time to check for charge
% leakage
% 
% --> Read Cycles vs. Resistance, evaluates how repeated reading affects
% retention

% * Power consumption (during swithcing)
%
% * Switching speed
%
% --> R - t plot, tracks resistance change
%
% --> I - t plot, measures transition time

%% B-V electrochemical equation

% Find/Plot the charge density v.s. position
dfplot.rhoxFxVx(sol);

% plot umax/uo v.s. ions density
% plot umax/uo v.s. thickness

% umax/u0 = (u x Emax)

%% - - - - - - - - - - Anion/Cation Equilibrate Density (离子平衡浓度)
dfplot.acx(soleq.ion);

%% - - - - - - - - - - Boundary voltage (E)
dfplot.Vx(soleq.ion);

%% - - - - - - - - - - Energy level and band
dfplot.ELx(soleq.ion);

%% - - - - - - - - - - Energy level, band and charge carriers density
dfplot.ELxnpxacx(soleq.ion);

%% - - - - - - - - - - Test

% time v.s. ions density boundary

% time v.s. overpotential (butler-volmer voltage, E - Eeq)

%% - - - - - - - - - - time v.s. total ionic charge
dfplot.Qt(soleq.ion, 900, 920);
