input = 'Input_files/pmpi_only.csv';
par_a = pc(input);
par_a.tmesh_type = 'linear';
par_a.j0 = 0;
par_a = refresh_device(par_a);
soleq_a = equilibrate(par_a);
sol_a = doCV(soleq_a.ion, 0, 0, -1, 1, 1e-1, 1, 500);

%% - - - - - - - - - -

J = dfana.calcJ(sol_a);
Vapp = dfana.calcVapp(sol_a);
t = sol_a.t;
xmesh = sol_a.x;
xpos = 0;
ppos = getpointpos(xpos, xmesh);

dfplot_ionic.c2c(J, Vapp, ppos, 1);
