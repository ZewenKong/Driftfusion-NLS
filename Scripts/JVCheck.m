% par = pc('Input_files/EnergyOffsetSweepParameters_v5_undoped_SAM.csv');
parr = pc('Input_files/1_layer_single_carrier.csv');
eqm = equilibrate(parr);
CV_sol_ion = doCV(eqm.ion, 1, 0, 1, -1, 1e-4, 1, 401);
% CV_sol_el = doCV(eqm.el, 1, -0.3, 1.7, -0.3, 1e-4, 1, 401);
%%
Vapp = dfana.calcVapp(CV_sol_ion);
Jion = dfana.calcJ(CV_sol_ion);

figure(666)
% plot(Vapp, Jion.tot(:,1));
plot(Vapp, Jion.tot(:,1))