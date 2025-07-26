classdef df_analysis

    % DRIFTFUSION analysis class- contains multiple methods for calculating
    % outputs using the solution obtained from DF.

    % LICENSE
    % Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
    % Imperial College London
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU Affero General Public License as published
    % by the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.

    % ADD-ONs: Modelling - ionic charge transport (validation)
    % func 'maxFatIF()', max potential at interface
    % func 'normdmu()', normalised ionic mobility

    % ADD-ONs: Modelling - Butler-Volmer calculation
    % func 'eqdens()', equilibrium ion density
    % func 'denstochemact()'. charge carroer density (number) to chemical activity
    % func 'nernst()', nernst equation
    % func 'butlervolmer()', butler-volmer equation

    % ADD-ONs: Calculated data process (save)
    % func 'process_calcdata()', process calculated data

    % ADD-ONs: special time point analysis (max J)
    % func 'maxJAndCorr()', find max J value and its corresponding value (same voltage)
    % func 'tPntOfSpJs()', find time points of these special value (based on J)
    % func 'spPntTosoleq()', use the speical points as the equilibrium function (max and its corresponding)

    % ADD-ONs: find the time point of special J or special V
    % func 'tOfSpJ()'
    % func 'tOfSpV()'

    % ADD-ONs: the pure ionic charge density (number)
    % func 'ionicrhocalc()'

    % - - - - - - - - - - CODE START - - - - - - - - - -

    methods (Static)

        function [u, t, x, par, dev, n, p, a, c, V] = splitsol(sol)
            % splits solution into useful outputs
            u = sol.u;
            t = sol.t(1:size(u, 1));
            x = sol.x;
            par = sol.par;
            dev = par.dev;
            % split the solution into its component parts (e.g. electrons, holes and efield)
            V = u(:, :, 1);
            n = u(:, :, 2);
            p = u(:, :, 3);

            switch par.N_ionic_species
                case 0
                    c = zeros(length(t), length(x));
                    a = zeros(length(t), length(x));
                    dev.Ncat = zeros(1, length(x));
                    dev.Nani = zeros(1, length(x));
                    par.dev.Ncat = zeros(1, length(x));
                    par.dev.Nani = zeros(1, length(x));
                    par.dev_sub.Ncat = zeros(1, length(x) - 1);
                    par.dev_sub.Nani = zeros(1, length(x) - 1);
                case 1
                    c = u(:, :, 4);
                    a = zeros(length(t), length(x));
                    dev.Nani = zeros(1, length(x));
                    par.dev.Nani = zeros(1, length(x));
                    par.dev_sub.Nani = zeros(1, length(x) - 1);
                case 2
                    c = u(:, :, 4);
                    a = u(:, :, 5);
            end

        end

        function [Ecb, Evb, Efn, Efp] = calcEnergies(sol)

            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            % u is the solution structure
            Ecb = dev.Phi_EA - V; % conduction band potential
            Evb = dev.Phi_IP - V; % valence band potential
            Efn = zeros(size(n, 1), size(n, 2));
            Efp = zeros(size(n, 1), size(n, 2));

            switch par.prob_distro_function
                case 'Fermi'

                    for i = 1:size(n, 1) % time

                        for j = 1:size(n, 2) % position
                            Efn(i, j) = distro_fun.Efn_fd_fun(n(i, j), dev.Efn(j, :), dev.n_fd(j, :));
                            Efp(i, j) = distro_fun.Efp_fd_fun(p(i, j), dev.Efp(j, :), dev.p_fd(j, :));
                        end

                    end

                    Efn = Efn - V;
                    Efp = Efp - V;

                case 'Blakemore'
                    Efn = real(Ecb + (par.kB * par.T / par.q) * log(n ./ (dev.Nc - par.gamma * n)));
                    Efp = real(Evb - (par.kB * par.T / par.q) * log(p ./ (dev.Nv - par.gamma * p)));

                case 'Boltz'
                    Efn = real(Ecb + (par.kB * par.T / par.q) * log(n ./ dev.Nc)); % electron quasi-Fermi level
                    Efp = real(Evb - (par.kB * par.T / par.q) * log(p ./ dev.Nv)); % hole quasi-Fermi level
            end

        end

        function [J, j, x] = calcJ(sol)
            % current, J and flux, j calculation from continuity equations
            % obtain SOL components for easy referencing

            [u, t, xmesh, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);

            n_sub = getvar_sub(n);
            p_sub = getvar_sub(p);
            a_sub = getvar_sub(a);
            c_sub = getvar_sub(c);

            x = par.x_sub;
            [~, ~, g] = df_analysis.calcg(sol);

            [~, dndt] = gradient(n_sub, x, t);
            [~, dpdt] = gradient(p_sub, x, t);
            [~, dadt] = gradient(a_sub, x, t);
            [~, dcdt] = gradient(c_sub, x, t);

            % recombination
            r = df_analysis.calcr(sol, "sub");

            djndx = -dndt + g - r.tot;
            djpdx = -dpdt + g - r.tot;
            djadx = -dadt; % add source terms as necessary
            djcdx = -dcdt; % add source terms as necessary

            deltajn = cumtrapz(x, djndx, 2);
            deltajp = cumtrapz(x, djpdx, 2);
            deltaja = cumtrapz(x, djadx, 2);
            deltajc = cumtrapz(x, djcdx, 2);

            % currents from the boundaries
            jn_l = -par.sn_l * (n(:, 1) - par.n0_l);
            jn_r = par.sn_r * (n(:, end) - par.n0_r);

            jp_l = -par.sp_l * (p(:, 1) - par.p0_l);
            jp_r = par.sp_r * (p(:, end) - par.p0_r);

            jc_l = 0;
            jc_r = 0;

            ja_l = 0;
            ja_r = 0;

            % calculate total electron and hole currents from fluxes
            % use the minority carrier flux as the boundary condition
            if par.p0_l == par.n0_l && par.n0_r == par.p0_r
                % intrinsic both sides then integrate from side with lowest density
                if par.n0_r > par.n0_l
                    j.n = jn_l + deltajn;
                else
                    j.n = jn_r + (deltajn - deltajn(:, end));
                end

                if par.p0_r > par.p0_l
                    j.p = jp_l + deltajp;
                else
                    j.p = jp_r + (deltajp - deltajp(:, end));
                end

            elseif par.p0_l >= par.n0_l && par.n0_r >= par.p0_r
                % p-type left boundary, n-type right boundary
                j.n = jn_l + deltajn;
                j.p = jp_r + (deltajp - deltajp(:, end));
            elseif par.n0_l >= par.n0_r && par.p0_r >= par.n0_r
                % n-type left boundary, p-type right boundary
                j.n = jn_r + (deltajn - deltajn(:, end));
                j.p = jp_l + deltajp;
            elseif par.p0_l >= par.n0_l && par.p0_r >= par.n0_r ...
                    || par.n0_l >= par.p0_l && par.n0_r >= par.p0_r
                % p-type both boundaries or n-type both boundaries
                j.n = jn_l + deltajn;
                j.p = jp_l + deltajp;
            end

            j.c = jc_l + deltajc;
            j.a = ja_l + deltaja;

            % apply switches and accelerators
            j.n = par.mobset * j.n;
            j.p = par.mobset * j.p;
            j.c = par.mobseti * par.K_c * j.c;
            j.a = par.mobseti * par.K_a * j.a;

            % displacement flux
            FV_sub = df_analysis.calcF(sol, "sub");

            [~, FV_sub_dt] = gradient(FV_sub, x, t);
            j.disp = par.epp0 .* par.dev_sub.epp .* FV_sub_dt;

            J.n = j.n * -par.e;
            J.p = j.p * par.e;
            J.c = j.c * par.z_c * par.e;
            J.a = j.a * par.z_a * par.e;
            J.disp = j.disp * par.e;

            % total current
            J.tot = J.n + J.p + J.a + J.c + J.disp;
        end

        function [g1, g2, g] = calcg(sol)
            [~, tmesh, ~, par, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);

            % generation function
            switch par.g1_fun_type
                case 'constant'
                    g1 = repmat(par.int1 .* par.gx1, length(tmesh), 1);
                otherwise
                    g1_fun = fun_gen(par.g1_fun_type);
                    g1 = g1_fun(par.g1_fun_arg, tmesh') * par.gx1;
            end

            switch par.g2_fun_type
                case 'constant'
                    g2 = repmat(par.int2 .* par.gx2, length(tmesh), 1);
                otherwise
                    g2_fun = fun_gen(par.g2_fun_type);
                    g2 = g2_fun(par.g2_fun_arg, tmesh') * par.gx2;
            end

            g = g1 + g2;
        end

        function Pin = calcPin(sol)
            % incident optical power density
            % note this integrates across the available spectrum within AM15.xls
            if strcmp(sol.par.optical_model, 'Beer-Lambert')
                AM15_data = readtable('AM15.xls', 'VariableNamingRule', 'preserve');
                Pin = 1e-3 * trapz(AM15_data.(1), AM15_data.(2));
            else
                warning('No incident photon spectrum available, assuming Pin = 0.1 W cm-2')
                Pin = 0.1;
            end

        end

        % calculate the recombination rate on i-half mesh, obtain SOL components for easy referencing
        function [r, ns, ps, alpha_xn, beta_xp] = calcr(sol, mesh_option)
            % MESH_OPTION = "whole" for input mesh or "sub" for subinetrval mesh
            [u, t, x_input, par, ~, n, p, a, c, V] = df_analysis.splitsol(sol);

            switch mesh_option
                case "whole"
                    dev = par.dev;
                    x = x_input;
                    n = u(:, :, 2);
                    p = u(:, :, 3);
                case "sub"
                    dev = par.dev_sub;
                    x = par.x_sub;
                    n = getvar_sub(u(:, :, 2));
                    p = getvar_sub(u(:, :, 3));
            end

            dVdx = zeros(length(t), length(x));

            for i = 1:length(t)
                [~, dVdx(i, :)] = pdeval(0, x_input, V(i, :), x);
            end

            vsr_zone = repmat(dev.vsr_zone, length(t), 1);
            srh_zone = repmat(dev.srh_zone, length(t), 1);

            xprime_n = dev.xprime_n;
            xprime_p = dev.xprime_p;
            sign_xn = repmat(dev.sign_xn, length(t), 1); % 1 if xn increasing, -1 if decreasing wrt x
            sign_xp = repmat(dev.sign_xp, length(t), 1); % 1 if xp increasing, -1 if decreasing wrt x
            alpha0_xn = repmat(dev.alpha0_xn, length(t), 1);
            beta0_xp = repmat(dev.beta0_xp, length(t), 1);

            alpha_xn = (sign_xn .* par.q .* dVdx ./ (par.kB * par.T)) + alpha0_xn;
            beta_xp = (sign_xp .* par.q .* -dVdx ./ (par.kB * par.T)) + beta0_xp;

            % band-to-band
            r.btb = dev.B .* (n .* p - dev.ni .^ 2);

            % bulk SRH
            r.srh = srh_zone .* (n .* p - dev.ni .^ 2) ...
                ./ (dev.taun .* (p + dev.pt) + dev.taup .* (n + dev.nt));

            % volumetric surface SRH
            ns = n .* exp(-alpha_xn .* xprime_n); % projected electron surface density
            ps = p .* exp(-beta_xp .* xprime_p); % projected hole surface density
            r.vsr = vsr_zone .* (ns .* ps - dev.nt .* dev.pt) ...
                ./ (dev.taun_vsr .* (ps + dev.pt) + dev.taup_vsr .* (ns + dev.nt));
            % system boundary surface recombination i.e. minority carrier currents
            r.tot = r.btb + r.srh + r.vsr; % total

        end

        % calculates the absolute surface recombination flux for system boundaries
        function j_surf_rec = calcj_surf_rec(sol)

            [u, t, x_input, par, ~, n, p, a, c, V] = df_analysis.splitsol(sol);

            % absolute fluxes at the boundaries
            [~, j, ~] = df_analysis.calcJ(sol);
            jn_l = abs(j.n(:, 1));
            jn_r = abs(j.n(:, end));

            jp_l = abs(j.p(:, 1));
            jp_r = abs(j.p(:, end));

            j_surf_rec.n_l = zeros(1, length(t));
            j_surf_rec.p_l = zeros(1, length(t));
            j_surf_rec.n_r = zeros(1, length(t));
            j_surf_rec.p_r = zeros(1, length(t));

            if par.p0_l == par.n0_l && par.n0_r == par.p0_r
                % intrinsic both sides then can be either?
                j_surf_rec.l = jn_l;
                j_surf_rec.r = jn_r;

            elseif par.p0_l >= par.n0_l && par.n0_r >= par.p0_r
                % p-type left boundary, n-type right boundary
                j_surf_rec.l = jn_l;
                j_surf_rec.r = jp_r;

            elseif par.n0_l >= par.n0_r && par.p0_r >= par.n0_r
                % n-type left boundary, p-type right boundary
                j_surf_rec.l = jp_l;
                j_surf_rec.r = jn_r;

            elseif par.p0_l >= par.n0_l && par.p0_r >= par.n0_r
                j_surf_rec.l = jn_l;
                j_surf_rec.r = jn_r;

            elseif par.n0_l >= par.p0_l && par.n0_r >= par.p0_r
                % p-type both boundaries or n-type both boundaries
                j_surf_rec.l = jp_l;
                j_surf_rec.r = jp_r;
            end

            j_surf_rec.tot = j_surf_rec.l + j_surf_rec.r;
        end

        function j_surf_rec = calcj_surf_rec_lucy(sol)
            [u, t, x_input, par, ~, n, p, a, c, V] = df_analysis.splitsol(sol);
            [~, j, ~] = df_analysis.calcJ(sol);

            jn_l = j.n(:, 1);
            jn_r = j.n(:, end);

            jp_l = j.p(:, 1);
            jp_r = j.p(:, end);

            for i = 1:length(t)
                j_surf_rec.l(i) = jn_l(i);
                j_surf_rec.r(i) = jp_r(i);
                j_surf_rec.tot(i) = abs(j_surf_rec.l(i)) + abs(j_surf_rec.r(i)); %+2*j_gen(1);
            end

        end

        % calculates drift and diffusion currents at every point and all times -
        function [Jdd, jdd, xout] = calcJdd(sol)
            % NOTE: UNRELIABLE FOR TOTAL CURRENT as errors in the calculation of the
            % spatial gradients mean that the currents do not cancel properly
            % obtain SOL components for easy referencing
            [~, t, x, par, ~, n, p, a, c, V] = df_analysis.splitsol(sol);
            xout = par.x_sub;
            dev = par.dev_sub;

            % property matrices
            eppmat = dev.epp;
            mu_n_mat = dev.mu_n;
            mu_p_mat = dev.mu_p;
            mu_cat = dev.mu_c;
            mu_ani = dev.mu_a;
            gradEA_mat = dev.gradEA;
            gradIP_mat = dev.gradIP;
            gradNc_mat = dev.gradNc;
            gradNv_mat = dev.gradNv;
            Nc_mat = dev.Nc;
            Nv_mat = dev.Nv;

            V_sub = zeros(length(t), length(xout));
            n_sub = zeros(length(t), length(xout));
            p_sub = zeros(length(t), length(xout));
            a_sub = zeros(length(t), length(xout));
            c_sub = zeros(length(t), length(xout));

            dVdx = zeros(length(t), length(xout));
            dndx = zeros(length(t), length(xout));
            dpdx = zeros(length(t), length(xout));
            dadx = zeros(length(t), length(xout));
            dcdx = zeros(length(t), length(xout));

            % avoid PDEVAL for faster calculation
            % obtain variables and gradients on sub-interval mesh
            for i = 1:length(t)
                V_sub(i, :) = 0.5 * (V(i, 2:end) + V(i, 1:end - 1));
                n_sub(i, :) = 0.5 * (n(i, 2:end) + n(i, 1:end - 1));
                p_sub(i, :) = 0.5 * (p(i, 2:end) + p(i, 1:end - 1));
                c_sub(i, :) = 0.5 * (c(i, 2:end) + c(i, 1:end - 1));
                a_sub(i, :) = 0.5 * (a(i, 2:end) + a(i, 1:end - 1));

                dVdx(i, :) = (V(i, 2:end) - V(i, 1:end - 1)) ./ (x(2:end) - x(1:end - 1));
                dndx(i, :) = (n(i, 2:end) - n(i, 1:end - 1)) ./ (x(2:end) - x(1:end - 1));
                dpdx(i, :) = (p(i, 2:end) - p(i, 1:end - 1)) ./ (x(2:end) - x(1:end - 1));
                dcdx(i, :) = (c(i, 2:end) - c(i, 1:end - 1)) ./ (x(2:end) - x(1:end - 1));
                dadx(i, :) = (a(i, 2:end) - a(i, 1:end - 1)) ./ (x(2:end) - x(1:end - 1));
            end

            % diffusion coefficients
            switch par.prob_distro_function
                case 'Fermi'

                    for jj = 1:length(x)
                        Dn_mat(i, jj) = distro_fun.D(n(i, jj), dev.Dnfun(jj, :), dev.n_fd(jj, :));
                        Dp_mat(i, jj) = distro_fun.D(p(i, jj), dev.Dpfun(jj, :), dev.p_fd(jj, :));
                    end

                case 'Blakemore'
                    Dn_mat = mu_n_mat .* par.kB .* par.T .* (Nc_mat ./ (Nc_mat - par.gamma .* n_sub));
                    Dp_mat = mu_p_mat .* par.kB .* par.T .* (Nv_mat ./ (Nv_mat - par.gamma .* p_sub));

                case 'Boltz'
                    Dn_mat = mu_n_mat * par.kB * par.T;
                    Dp_mat = mu_p_mat * par.kB * par.T;
            end

            % particle fluxes (remember F = -dVdx)
            jdd.ndrift = mu_n_mat .* n_sub .* (dVdx - gradEA_mat);
            jdd.ndiff = -Dn_mat .* (dndx - ((n_sub ./ Nc_mat) .* gradNc_mat));
            jdd.pdrift = mu_p_mat .* p_sub .* (-dVdx + gradIP_mat);
            jdd.pdiff = -Dp_mat .* (dpdx - ((p_sub ./ Nv_mat) .* gradNv_mat));

            switch par.N_ionic_species
                case 0
                    jdd.cdrift = zeros(length(t), length(xout));
                    jdd.cdiff = zeros(length(t), length(xout));
                    jdd.adrift = zeros(length(t), length(xout));
                    jdd.adiff = zeros(length(t), length(xout));
                case 1
                    jdd.cdrift = mu_cat .* c_sub .* -dVdx;
                    jdd.cdiff = -mu_cat .* par.kB * par.T .* dcdx;
                    jdd.adrift = zeros(length(t), length(xout));
                    jdd.adiff = zeros(length(t), length(xout));
                case 2
                    jdd.cdrift = mu_cat .* c_sub .* -dVdx;
                    jdd.cdiff = -mu_cat .* par.kB * par.T .* dcdx;
                    jdd.adrift = -mu_ani .* a_sub .* -dVdx;
                    jdd.adiff = -mu_ani .* par.kB * par.T .* dadx;
            end

            % note these have a negative sign compared with the expressions in dfNLSDE
            jdd.n = par.mobset * (jdd.ndrift + jdd.ndiff);
            jdd.p = par.mobset * (jdd.pdrift + jdd.pdiff);
            jdd.a = par.mobseti * (jdd.adrift + jdd.adiff);
            jdd.c = par.mobseti * (jdd.cdrift + jdd.cdiff);

            % displacement current
            [~, dFdt] = gradient(-dVdx, xout, t);
            j.disp = par.epp0 .* eppmat .* dFdt;

            jdd.disp = j.disp;
            % the total flux here includes the sign of the carrier
            jdd.tot = -jdd.n + jdd.p - jdd.a + jdd.c + jdd.disp;

            Jdd.ndrift = jdd.ndrift * par.e;
            Jdd.ndiff = jdd.ndiff * par.e;
            Jdd.pdrift = jdd.pdrift * par.e;
            Jdd.pdiff = jdd.pdiff * par.e;
            Jdd.adrift = jdd.adrift * par.e;
            Jdd.adiff = jdd.adiff * par.e;
            Jdd.cdrift = jdd.cdrift * par.e;
            Jdd.cdiff = jdd.cdiff * par.e;

            Jdd.n = par.e * jdd.n;
            Jdd.p = par.e * jdd.p;
            Jdd.a = par.e * jdd.a;
            Jdd.c = par.e * jdd.c;
            Jdd.disp = par.e * jdd.disp;
            Jdd.tot = par.e * jdd.tot;
        end

        function [FV, Frho] = calcF(sol, mesh_option)
            % dlectric field caculation
            %
            % FV = Field calculated from the gradient of the potential
            % Frho = Field calculated from integrated space charge density

            [u, t, x_whole, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);

            switch mesh_option
                case "whole"
                    x = x_whole;
                case "sub"
                    x = par.x_sub;
            end

            for i = 1:length(t)
                [~, dVdx(i, :)] = pdeval(0, x_whole, V(i, :), x);
            end

            FV = -dVdx;

            if nargout > 1
                rho = df_analysis.calcrho(sol, "sub");
                Frho = cumtrapz(x, rho, 2) ./ (par.dev_sub.epp .* par.epp0) + FV(:, 1);
            end

        end

        function rho = calcrho(sol, mesh_option)

            % calculates the space charge density

            [u, t, x, par, dev_in, n_whole, p_whole, a_whole, c_whole, V_whole] = df_analysis.splitsol(sol);

            switch mesh_option
                case "whole"
                    dev = par.dev;
                    n = n_whole;
                    p = p_whole;
                    a = a_whole;
                    c = c_whole;
                case "sub"
                    dev = par.dev_sub;
                    n = getvar_sub(n_whole);
                    p = getvar_sub(p_whole);
                    a = getvar_sub(a_whole);
                    c = getvar_sub(c_whole);
            end

            NA = repmat(dev.NA, length(t), 1);
            ND = repmat(dev.ND, length(t), 1);
            Nani = repmat(dev.Nani, length(t), 1);
            Ncat = repmat(dev.Ncat, length(t), 1);

            % charge density
            rho = -n + p - NA + ND + par.z_a * a + par.z_c * c - par.z_a * Nani - par.z_c * Ncat;
        end

        function Vapp = calcVapp(sol)
            [~, t, ~, par, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);

            switch par.V_fun_type
                case 'constant'
                    Vapp = ones(1, length(t)) * par.V_fun_arg(1);
                otherwise
                    Vapp_fun = fun_gen(par.V_fun_type);
                    Vapp = Vapp_fun(par.V_fun_arg, t);
            end

        end

        function stats = JVstats(JVsol)

            % a function to pull statistics from a JV sweep using DOJV
            % JVsol - a solution from DOJV

            if isfield(JVsol, 'ill')

                if isfield(JVsol.ill, 'f')
                    Vapp = df_analysis.calcVapp(JVsol.ill.f);
                    Vapp = Vapp';
                    J = df_analysis.calcJ(JVsol.ill.f);

                    try
                        stats.Jsc_f = interp1(Vapp, J.tot(:, end), 0);
                    catch
                        warning('No Jsc available- Vapp must pass through 0 to obtain Jsc')
                        stats.Jsc_f = 0;
                    end

                    try
                        stats.Voc_f = interp1(J.tot(:, end), Vapp, 0);
                    catch
                        warning('No Voc available- try increasing applied voltage range')
                        stats.Voc_f = 0;
                    end

                    if stats.Jsc_f ~= 0 && stats.Voc_f ~= 0
                        pow_f = J.tot(:, end) .* Vapp;
                        stats.mpp_f = min(pow_f);
                        stats.mppV_f = Vapp(pow_f == stats.mpp_f);
                        stats.FF_f = stats.mpp_f / (stats.Jsc_f * stats.Voc_f);
                    end

                    % hysteresis index
                    A_f = abs(trapz(Vapp(Vapp >= 0 & Vapp <= stats.Voc_f), J.tot(Vapp >= 0 & Vapp <= stats.Voc_f, end)));

                else
                    stats.Jsc_f = nan;
                    stats.Voc_f = nan;
                    stats.mpp_f = nan;
                    stats.FF_f = nan;
                end

                if isfield(JVsol.ill, 'r')
                    Vapp = df_analysis.calcVapp(JVsol.ill.r);
                    Vapp = Vapp';
                    J = df_analysis.calcJ(JVsol.ill.r);

                    try
                        stats.Jsc_r = interp1(Vapp, J.tot(:, end), 0);
                    catch
                        warning('No Jsc available- Vapp must pass through 0 to obtain Jsc')
                        stats.Jsc_r = 0;
                    end

                    try
                        stats.Voc_r = interp1(J.tot(:, end), Vapp, 0);
                    catch
                        warning('No Voc available- try increasing applied voltage range')
                        stats.Voc_r = 0;
                    end

                    if stats.Jsc_r ~= 0 && stats.Voc_r ~= 0
                        pow_r = J.tot(:, end) .* Vapp;
                        stats.mpp_r = min(pow_r);
                        stats.mppV_r = Vapp(pow_r == stats.mpp_r);
                        stats.FF_r = stats.mpp_r / (stats.Jsc_r * stats.Voc_r);
                    end

                    % hysteresis factor
                    A_r = abs(trapz(Vapp(Vapp >= 0 & Vapp <= stats.Voc_r), J.tot(Vapp >= 0 & Vapp <= stats.Voc_r, end)));

                    % sign to identify inverted hysteresis
                    if A_r >= A_f
                        B = 1;
                    elseif A_r < A_f
                        B = -1;
                    end

                    stats.HF = B * abs((A_r - A_f) / A_r);

                else
                    stats.Jsc_r = NaN;
                    stats.Voc_r = NaN;
                    stats.mpp_r = NaN;
                    stats.FF_r = NaN;
                    stats.HF = NaN;
                end

            else
            end

        end

        function value = calcPLt(sol)
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);

            Bmat = dev.B;
            value = trapz(x, (dev.B .* (n .* p - dev.ni .^ 2)), 2);
        end

        function DeltaQFL = calcDeltaQFL(sol)
            % get QFLs
            [~, ~, Efn, Efp] = df_analysis.calcEnergies(sol);
            par = sol.par;

            if par.p0_l >= par.n0_l && par.n0_r >= par.p0_r
                % p-type left boundary, n-type right boundary
                DeltaQFL = Efn(:, end) - Efp(:, 1);

            elseif par.n0_l >= par.n0_r && par.p0_r >= par.n0_r
                % n-type left boundary, p-type right boundary
                DeltaQFL = Efn(:, 1) - Efp(:, end);

            else
                % if all equal the choose arbitrary boundaries
                DeltaQFL = Efn(:, end) - Efp(:, 1);
            end

        end

        function deltaV = deltaVt(sol, p1, p2)
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            % calculates the electrostatic potential difference as a function of time between two points P1 and P2
            deltaV = V(:, p1) - V(:, p2);
        end

        function sigma = calcsigma(sol)
            % calculates the integrated space charge density
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            rho = df_analysis.calcrho(sol, "whole");
            sigma = trapz(x, rho, 2);
        end

        function sigma_ion = calcsigma_ion(sol)
            % calculates the integrated space charge density
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            rho_ion = c - a;
            sigma_ion = trapz(x, rho_ion, 2);
        end

        function Fion = calcFion(sol)
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);

            rhoion = c - dev.Ncat - a + dev.Nani;
            Fion = cumtrapz(x, rhoion, 2) ./ (dev.epp * par.epp0);

        end

        function Vion = calcVion(sol)
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);

            Fion = df_analysis.calcFion(sol);
            Vion = -cumtrapz(x, Fion, 2);
        end

        function [U, Ux] = pdentrp(singular, m, xL, uL, xR, uR, xout)

            % PDENTRP  Interpolation helper function for PDEPE.
            %   [U,UX] = PDENTRP(M,XL,UL,XR,UR,XOUT) uses solution values UL at XL and UR at XR
            %   for successive mesh points XL < XR to interpolate the solution values U and
            %   the partial derivative with respect to x, UX, at arguments XOUT(i) with
            %   XL <= XOUT(i) <= XR.  UL and UR are column vectors. Column i of the output
            %   arrays U, UX correspond to XOUT(i).
            %
            %   See also PDEPE, PDEVAL, PDEODES.
            %
            %   Lawrence F. Shampine and Jacek Kierzenka
            %   Copyright 1984-2013 The MathWorks, Inc.
            %   $Revision: 1.5.4.4.54.1 $  $Date: 2013/09/27 03:10:22 $

            xout = xout(:)';
            nout = length(xout);

            U = uL(:, ones(1, nout));
            Ux = zeros(size(U));

            uRL = uR - uL;

            % use singular interpolant on all subintervals
            if singular
                U = U + uRL * ((xout .^ 2 - xL ^ 2) / (xR ^ 2 - xL ^ 2));
                Ux = uRL * (2 * xout / (xR ^ 2 - xL ^ 2));
            else

                switch m
                    case 0
                        U = U + uRL * ((xout - xL) / (xR - xL));
                        Ux = uRL * (ones(1, nout) / (xR - xL));
                    case 1
                        U = U + uRL * (log(xout / xL) / log(xR / xL));
                        Ux = uRL * ((1 ./ xout) / log(xR / xL));
                    case 2
                        U = U + uRL * ((xR ./ xout) .* ((xout - xL) / (xR - xL)));
                        Ux = uRL * ((xR ./ xout) .* (xL ./ xout) / (xR - xL));
                end

            end

        end

        % ADD-ONs: Modelling - ionic charge transport

        % find the max electrical field value (point) at interface
        function [F_idx, maxF_idx] = maxFatIF(sol)
            [u, t, x_whole, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            F = df_analysis.calcF(sol, "whole"); % the potential through the device
            tot_points = cumsum(par.layer_points); % accumulate the position points
            boundary_idx = tot_points(length(par.layer_points) - 2); % points at the boundary (interface) index
            F_idx = F(:, boundary_idx); % all potential (over time) at that boundary point
            maxF_idx = max(F_idx); % max potential at that point
        end

        % calculate normalised ionic mobility
        function mun = normdmu(Fmax)
            kB = 1.38e-23; % boltzmann const.
            T = 298; % room temperature
            d = 6.32e-8; % 2DRP lattice const. (cm)

            % https://pubs.aip.org/aip/jap/article-abstract/128/17/175101/1063040/
            % Rashba-band-splitting-in-two-dimensional?redirectedFrom=fulltext
            % https://advanced.onlinelibrary.wiley.com/doi/10.1002/adfm.201903293

            z = 1; % charge number of the ion
            e = 1.6e-19; % elementary charge
            F = abs(Fmax); % the max energy for each sample

            % probability term
            Ea = 0.45 * e; % activation energy (eV), approximate value
            v0 = 5e12; % phonon frequency (General value)
            c = 1; % for a 1D transport, approximately equal to unity

            % equation
            P = c * v0 * exp(- (Ea / (kB * T))); % transition probability
            vD = d * P * (exp((z * e * d * F) / (2 * kB * T)) ...
                - exp(- (z * e * d * F) / (2 * kB * T)));
            vD0 = d * P * ((z * e * d * F) / (kB * T));
            mu = vD ./ F; % mobility (drfit velocity divide field potential)
            mu0 = vD0 ./ F; % mobility at low electric field (dF/kBT << 1)
            mun = mu ./ mu0;
        end

        % ADD-ONs: Butler-Volmer calculation

        % equilibrium ionic charge carrier density calculation (anion here)
        function a_eq = eqdens(E_b, E_st, R, T, F, A, z)
            % equilibrium anion (Ag+) density calculation
            a_eq = (exp((E_b - E_st) / (- (R * T) / (F * z))) * A) / (1e3);
        end

        % charge carrier density to chemical activity
        function a = denstochemact(density, A)
            a = (density * 1e3) / A;
            % for a given dissolved species, its chemical activity (a)
            % is the product of its activity coefficient (γ) by its molar (mol/L solution),
            % or molal (mol/kg water), concentration (C): a = γ C
        end

        % Nernst equation
        function E_eq = nernst(a, E_st, R, T, F, z)
            E_eq = E_st - ((R * T) / (F * z)) * log(a);
            % a, the chemical activity
            % E_st, the standard potential (Ag/AgI = 0.152)
            % R, Universal ideal gas constant
            % T, Temperature
            % F, Faraday constant
            % z, number of charge transfer (1)
        end

        % Butler-Volmer equation (ionic current flux)
        function j_bv = butlervolmer(j0, alpha_e, R, T, F, E, E_eq)

            j_temp = j0;
            R_series = 0.15; % series resistance 0.15 ohm/cm2
            eta_modified = (E - E_eq) - j_temp * R_series;
            eta_lim = 0.3; % eta limitation
            eta_soft = eta_modified ./ sqrt(1 + (eta_modified / eta_lim) .^ 2);
            j_bv = (- j0) * (exp((1 - alpha_e) * F * eta_soft / (R * T)) - exp((- alpha_e) * F * eta_soft / (R * T)));

            % j_bv = (- j0) * (exp((1 - alpha_e) * F * (E - E_eq) / (R * T)) - exp((- alpha_e) * F * (E - E_eq) / (R * T)));
            % j_bv = j0 * (exp((1 - alpha_e) * F * (E - E_eq) / (R * T)) - exp((- alpha_e) * F * (E - E_eq) / (R * T)));
            % j0, exchange current density (no bias)
            % * negative sign accounts for current direction being opposite to driftfusion direction sign
            % alpha_e, electrode (anodic/cathodic) charge transfer coefficient, normally 1/2
            % E is the electrode potential on the electrochemical scale
            % E_eq, is the equilibrium potential calculated by Nernst equation on the electrochemical scale
        end

        % ADD-ONs: calculated data process

        function process_calcdata(varargin)

            persistent calc; % persistent variable, the data value will be saved after function end (calculated results)

            if nargin == 1 && ischar(varargin{1}) && strcmp(varargin{1}, 'reset') % if input is "reset"
                calc = []; % reset the persistent variable
                disp('df_analysis.m - func process_calcdata: persistent calc has been reset');
                disp('-');
                return;
            end

            if isempty(calc) % if the calc = []
                calc = struct(); % initialise the structure
            end

            for i = 1:2:nargin % save the data accroding to the numerical field (字段)
                field_name = varargin{i}; % odd number-th field >> field name
                field_value = varargin{i + 1}; % even number-th field >> field value

                if ~isfield(calc, field_name)
                    calc.(field_name) = []; % initialise the field if it doesn't exist
                end

                calc.(field_name)(end + 1, 1) = field_value; % append the value to the corresponding field
            end

            assignin('base', 'calc', calc);
        end

        % ADD-ONs: special time point analysis (max J)

        % find max J and its corresponding value (same voltage)
        function [J_max, J_corr, idx_J_max, idx_J_corr] = maxJAndCorr(sol, xpos)
            J = df_analysis.calcJ(sol);
            Vapp = df_analysis.calcVapp(sol);
            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);
            J_temp = J.tot(:, ppos); % save as temp value
            Vapp_temp = Vapp;

            [J_max, idx_J_max] = max(J_temp); % get max J and its idx
            V_corr_J_max = Vapp_temp(idx_J_max); % get max J corresponding Vapp

            % filter the Vapp_temp, find the V values smaller than / equal to V_corr_J_max
            idx_filtered_arr = find(Vapp_temp <= V_corr_J_max);

            % transfer to the max diff array, get the max diff value
            [max_diff, idx_max_diff] = max(diff(idx_filtered_arr));

            % plus 1 to get the index of minuend, which is the index of corresponding J value
            idx_J_corr = idx_filtered_arr(idx_max_diff + 1);
            J_corr = J_temp(idx_J_corr);
        end

        % time points of special current fluxes
        function [t_J_max, t_J_corr] = tPntOfSpJs(sol, xpos)
            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);
            J = df_analysis.calcJ(sol);
            Vapp = df_analysis.calcVapp(sol);
            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);
            J_temp = J.tot(:, ppos);
            Vapp_temp = Vapp;

            [J_max, J_corr] = df_analysis.maxJAndCorr(sol, xpos);
            idx_J_max_1 = find(J_temp == J_max);
            idx_J_corr_1 = find(J_temp == J_corr);
            t_J_max = t(idx_J_max_1);
            t_J_corr = t(idx_J_corr_1);
        end

        % use the special / corresponding point as the equilibrium function
        function [selsoleq_1, selsoleq_2] = spPntTosoleq(sol, xpos) % selected (max J and its corr) point to equilibrium solution

            % solution processing
            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);
            J = df_analysis.calcJ(sol);
            Vapp = df_analysis.calcVapp(sol);

            % limitation
            [J_max, J_corr, idx_J_max, idx_J_corr] = df_analysis.maxJAndCorr(sol, xpos);
            [t_J_max, t_J_corr] = df_analysis.tPntOfSpJs(sol, xpos);
            idx_t_J_max = find(t == t_J_max);
            idx_t_J_corr = find(t == t_J_corr);

            % point of max J as new equilibrium solution
            selsoleq_1 = sol;
            selsoleq_1.u = sol.u(idx_t_J_max, :, :);
            selsoleq_1.t = sol.t(idx_t_J_max);
            V_corr_J_max = Vapp(idx_J_max);
            selsoleq_1.par.mobseti = 0; % * set mobseti = 0;
            selsoleq_1.par.j0 = 0;
            selsoleq_1.par = refresh_device(selsoleq_1.par);
            [~, sol_dwell_1] = ramped_step(selsoleq_1, V_corr_J_max, 0.1, 0.1); % ramp the solution
            sol_dwell_1.par.mobseti = 1; % set mobseti back
            sol_dwell_1.par = refresh_device(sol_dwell_1.par);
            selsoleq_1 = sol_dwell_1;

            % corrsponding point of max J as new equilibrium solution
            selsoleq_2 = sol;
            selsoleq_2.u = sol.u(idx_t_J_corr, :, :);
            selsoleq_2.t = sol.t(idx_t_J_corr);
            V_corr_J_corr = Vapp(idx_J_corr);
            selsoleq_2.par.mobseti = 0;
            selsoleq_2.par.j0 = 0;
            selsoleq_2.par = refresh_device(selsoleq_2.par);
            [~, sol_dwell_2] = ramped_step(selsoleq_2, V_corr_J_corr, 0.1, 0.1);
            sol_dwell_2.par.mobseti = 1;
            sol_dwell_2.par = refresh_device(sol_dwell_2.par);
            selsoleq_2 = sol_dwell_2;
        end

        % ADD-ONs: find the special point through J or V

        % find time point of special voltage point
        function t_sp = tOfSpV(sol, xpos, V_sp)
            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);
            Vapp = df_analysis.calcVapp(sol);
            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);
            V_temp = Vapp;
            % approximate 'find'
            [~, idx_V_sp] = min(abs(V_temp - V_sp)); % idx_V_sp = find(V_temp == V_sp);
            disp(['Sp V idx: ', num2str(idx_V_sp)]);
            t_sp = t(idx_V_sp);
        end

        % find time point of special J point
        function t_sp = tOfSpJ(sol, xpos, J_sp)
            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);
            J = df_analysis.calcJ(sol);
            J_temp = J.tot;
            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);
            [~, idx_J_sp] = min(abs(J_temp - J_sp));
            % disp(['Sp J idx: ', num2str(idx_J_sp)]);
            t_sp = t(idx_J_sp);
        end

        % ADD-ONs: find the pure ionic charge carriers density (number)

        function rho_ionic = rhoIonicCalc(sol, mesh_option)
            [u, t, x, par, dev_in, n_whole, p_whole, a_whole, c_whole, V_whole] = df_analysis.splitsol(sol);

            switch mesh_option
                case "whole"
                    dev = par.dev;
                    n = n_whole;
                    p = p_whole;
                    a = a_whole;
                    c = c_whole;
                case "sub"
                    dev = par.dev_sub;
                    n = getvar_sub(n_whole);
                    p = getvar_sub(p_whole);
                    a = getvar_sub(a_whole);
                    c = getvar_sub(c_whole);
            end

            NA = repmat(dev.NA, length(t), 1);
            ND = repmat(dev.ND, length(t), 1);
            Nani = repmat(dev.Nani, length(t), 1);
            Ncat = repmat(dev.Ncat, length(t), 1);
            rho_ionic = par.z_a * a + par.z_c * c - par.z_a * Nani - par.z_c * Ncat;
        end

        % calculate back
        % function eta_reverse = reversecalc(E, E_st, R, T, F, z, rho)
        %     eta_reverse = E - (E_st - ((R * T) / (F * z)) * log(rho));
        % end

    end

end
