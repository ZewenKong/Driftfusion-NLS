classdef df_plot
    % DRIFTFUSION Plotting class - contains methods for plotting
    %
    % List of available plots:

    % current related functions:
    % func 'JT()', currents as a function of time
    % ADD-ONs: func 'Jt_species()', sepecies' current v.s. time
    % func 'JV()', current - voltage curve using a solution from DOJV, now outdated - preferrable to use DOCV and df_plot.JTOTVAPP
    % func 'JVAPP()', current components as a function of the applied voltage at position defined by XPOS
    % func 'JTOTVAPP()', total current as a function of the applied voltage at position defined by XPOS
    % func 'LOGJVAPP()', current components as a function of the applied voltage using log y-axis at position defined by XPOS

    % charge density related functions:
    % func 'SIGMAT()', integrated charge density in cm-2 as a function of time
    % func 'QT()', charge density in Coulombs cm-2 integrated between within the range [X1, X2] as a function of time
    % ADD-ONs: func 'QionicT()', ionic charge density v.s. time
    % func 'QVAPP()', charge density in Coulombs cm-2 integrated between within the range [X1, X2] as a function of applied voltage
    % func 'RHOX(), volumetric charge density as a function of position
    % func 'DELTARHOX()', change in volumetric charge density as a function of position
    % func 'RHOXFXVX()', volumetric charge density, Electric field and Electrostatic potential as a function of position plot
    % func 'RHOXVX()', volumetric charge density and Electrostatic potential as a function of position plot

    % voltage related functions:
    % func 'VOCT()', open circuit voltage as a function of time
    % func 'VAPPT()', applied voltage as a function of time
    % func 'JVAPP()', current components as a function of the applied voltage at position defined by XPOS
    %
    % electrical field related functions:
    % func 'FT()', electric field as a function of time
    % func 'VXACX()', electrostatic potential and anion and cation densities, 2 panel
    % func 'VIONXACX()', electrostatic potential due to ionic charge and anion and cation densities, 2 panel
    % func 'FIONT()', electric field due to the ionic charge as a function of time

    % position related functions:
    % func 'XMESH()', plots the xmesh position vs point number
    % func 'JX()', total currents as a function of position
    % func 'jx()', carrier fluxes as a function of position
    % func 'JDDX()', drift and diffusion currents as a function of position
    % func 'VX()', electrostatic potential as a function of position
    % func 'NPX()', electron and hole densities as a function of position
    % func 'ACX()', anion and cation densities as a function of position
    % func 'GX()', generation rate as a function of position
    % func 'GXT()', generation rate as a function of position and time
    % func 'RX()', recombination rate components as a function of position

    % energy level related functions:
    % func 'ELX()', energy level diagram as a function of position
    % func 'ELXNPX()', energy level diagram, electron and hole densities
    % func 'ELXNPXACX()', energy level diagram, electron and hole densities, and anion and cation densities, 3 panel, materialed

    % others:
    % ADD-ONs: 'npxacx()', charge carriers species
    % func 'JRECVAPP()', recomonbination components as a function of the applied voltage
    % func 'PLT()', integrated radiative recombination rate as a function of time

    % ADD-ONs: func 'c2c()', cycle-to-cycle (C2C) characterisation
    % ADD-ONs: func 'CVappTrack()', current-Vapp track
    % ADD-ONs: func 'd2d_var()', device-to-device (D2D) characterisation with one variable
    % ADD-ONs: func 'd2d_vars()', D2D characterisation with two variables
    % ADD-ONs: func 'etaJ()', overpotential v.s. current (species), need pre-calculation

    % plotting functions that are a function of position can accept a time
    % array as the second argument- the procedure will loop and plot the
    % solution at multiple times.
    % The third optional argument defines the x-range.
    % For plotting functions that are a function of time, the second argument
    % is generally the position at which the value is taken- see the comments
    % of individual methods below for further details

    % LICENSE
    % Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
    % Imperial College London
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU Affero General Public License as published
    % by the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    %
    % - - - - - - - - - - CODE START - - - - - - - - - -

    methods (Static)

        % currents as a function of time
        function Jt(sol, xpos)
            % SOL = solution structure
            % XPOS = the readout position
            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);

            [J, j, xmesh] = df_analysis.calcJ(sol);
            ppos = getpointpos(xpos, xmesh);

            figure(2);
            plot(t, J.n(:, ppos), t, J.p(:, ppos), t, J.a(:, ppos), t, J.c(:, ppos), t, J.disp(:, ppos), t, J.tot(:, ppos));
            legend('Jn', 'Jp', 'Ja', 'Jc', 'Jdisp', 'Jtotal')
            xlabel('time [s]');
            ylabel('J [A cm^{-2}]');
            set(legend, 'FontSize', 16);
            set(legend, 'EdgeColor', [1 1 1]);
        end

        function Jtott(sol, xpos)
            % SOL = solution structure
            % XPOS = the readout position
            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);
            [J, j, xmesh] = df_analysis.calcJ(sol);
            ppos = getpointpos(xpos, xmesh);
            figure(2);
            plot(t, J.tot(:, ppos));
            xlabel('time [s]');
            ylabel('J [A cm^{-2}]');
        end

        % plots the current components
        function Jx(varargin)

            % VARARGIN = [SOL, TARR, XRANGE]
            % SOL = Solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [XMIN, XMAX]
            %
            % if no tarr, take the last point

            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            [J, j, x] = df_analysis.calcJ(sol);

            figure(3);
            df_plot.x2d(sol, x, {J.n, J.p, J.a, J.c, J.disp, J.tot}, ...
                {'Jn', 'Jp', 'Ja', 'Jc', 'Jdisp', 'Jtot'}, {'-', '-', '-', '-', '-', '-'}, ...
                'Current density [Acm-2]', tarr, xrange, 0, 0);
        end

        % plot the carrier fluxes
        function jx(varargin)

            % VARARGIN = [SOL, TARR, XRANGE]
            % SOL = Solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [XMIN, XMAX]

            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            [J, j, x] = df_analysis.calcJ(sol);
            figure(301);
            df_plot.x2d(sol, par.x_sub, {j.n, j.p, j.a, j.c, j.disp}, {'jn', 'jp', 'ja', 'jc', 'jdisp'}, ...
                {'-', '-', '-', '-', '-'}, 'Flux [cm-2 s-1]', tarr, xrange, 0, 0);
        end

        function JV(JV, option)
            % JV - a solution from doJV
            % OPTION - 1 = dark only, 2 = light only, 3 = dark & light
            % JV is a structure containing dark and illuminated JVs

            if option == 1 || option == 3
                J.dk.f = df_analysis.calcJ(JV.dk.f);
                Vapp.dk.f = df_analysis.calcVapp(JV.dk.f);
                J.dk.r = df_analysis.calcJ(JV.dk.r);
                Vapp.dk.r = df_analysis.calcVapp(JV.dk.r);

                figure(4)
                plot(Vapp.dk.f, J.dk.f.tot(:, end), '--', Vapp.dk.r, J.dk.r.tot(:, end));
                hold on
            end

            if option == 2 || option == 3

                J.ill.f = df_analysis.calcJ(JV.ill.f);
                Vapp.ill.f = df_analysis.calcVapp(JV.ill.f);
                J.ill.r = df_analysis.calcJ(JV.ill.r);
                Vapp.ill.r = df_analysis.calcVapp(JV.ill.r);

                figure(4)
                plot(Vapp.ill.f, J.ill.f.tot(:, end), '--') %,'Color', [0, 0.4470, 0.7410]);
                hold on
                plot(Vapp.ill.r, J.ill.r.tot(:, end)); %,'Color', [0, 0.4470, 0.7410]);
            end

            figure(4)
            xlabel('Applied voltage [V]')
            ylabel('Current density [Acm-2]');
            hold off
        end

        function Jddx(varargin)
            % drift and diffusion currents as a function of position
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            [Jdd, ~, x] = df_analysis.calcJdd(sol);

            figure(301);
            df_plot.x2d(sol, x, {Jdd.ndiff, Jdd.ndrift, Jdd.pdiff, Jdd.pdrift, ...
                                     Jdd.adiff, Jdd.adrift, Jdd.cdiff, Jdd.cdrift}, ...
                {'Jn,diff', 'Jn,drift', 'Jp,diff', 'Jp,drift', 'Ja,diff', 'Ja,drift', 'Jc,diff', 'Jc,drift'}, ...
                {'-', '-', '-', '-', '-', '-', '-', '-'}, 'Current density [Acm-2]', tarr, xrange, 0, 0);
        end

        function Voct(sol)
            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);
            Voc = df_analysis.calcDeltaQFL(sol);
            figure(6)
            plot(t, Voc)
            xlabel('Time [s]')
            ylabel('Voc [V]')
        end

        function PLt(sol)
            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);
            PL = df_analysis.calcPLt(sol);
            figure(7)
            plot(t, PL)
            xlabel('Time [s]')
            ylabel('PL [cm-2s-1]')
        end

        % difference in potential between the left and right boundary
        function Vappt(sol)
            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);
            Vapp = df_analysis.calcVapp(sol);

            figure(8)
            plot(t, Vapp);
            xlabel('Time [s]')
            ylabel('Vapp [V]')
        end

        % obtain point position from x position
        function JVapp(sol, xpos)

            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);
            J = df_analysis.calcJ(sol);
            Vapp = df_analysis.calcVapp(sol);

            figure(9)
            plot(Vapp, J.n(:, ppos), Vapp, J.p(:, ppos), Vapp, J.c(:, ppos), Vapp, J.a(:, ppos), Vapp, J.disp(:, ppos), Vapp, J.tot(:, ppos));
            legend('Jn', 'Jp', 'Jc', 'Ja', 'Jdisp', 'Jtotal')
            xlabel('Applied Voltage, Vapp [V]');
            ylabel('Current Density, J [A cm^{-2}]');
            set(legend, 'FontSize', 16);
            set(legend, 'EdgeColor', [1 1 1]);
        end

        % obtain point position from x position
        function JtotVapp(sol, xpos)

            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);
            J = df_analysis.calcJ(sol);
            Vapp = df_analysis.calcVapp(sol);

            figure(91)
            plot(Vapp, J.tot(:, ppos));
            xlabel('Applied Voltage, Vapp [V]');
            ylabel('Current Density, J [A cm^{-2}]');
            set(legend, 'FontSize', 16);
            set(legend, 'EdgeColor', [1 1 1]);
        end

        % obtain point position from x position
        function JtotVapp_mirror(sol, xpos)

            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);
            J = df_analysis.calcJ(sol);
            Vapp = df_analysis.calcVapp(sol);

            figure(91)
            plot(-Vapp, -J.tot(:, ppos));
            xlabel('Applied Voltage, Vapp [V]');
            ylabel('Current Density, J [A cm^{-2}]');
            set(legend, 'FontSize', 16);
            set(legend, 'EdgeColor', [1 1 1]);
        end

        % plot the log of the mod J
        function logJVapp(sol, xpos)

            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);
            J = df_analysis.calcJ(sol);
            Vapp = df_analysis.calcVapp(sol);

            figure(10)
            semilogy(Vapp, abs(J.tot(:, ppos)), Vapp, abs(J.n(:, ppos)), Vapp, abs(J.p(:, ppos)), Vapp, abs(J.a(:, ppos)), Vapp, abs(J.c(:, ppos)), Vapp, abs(J.disp(:, ppos)));
            xlabel('Vapp [V]');
            ylabel('|J| [A cm^{-2}]');
            legend('Jtot', 'Jn', 'Jp', 'Ja', 'Jc', 'Jdisp')
            set(legend, 'FontSize', 16);
            set(legend, 'EdgeColor', [1 1 1]);
        end

        function logJVapp3D(sol, xpos, ylogon)

            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);

            t = sol.t;
            J = df_analysis.calcJ(sol);
            Vapp = df_analysis.calcVapp(sol)';
            Jtot = J.tot(:, ppos);

            figure(11)
            surface('XData', [Vapp Vapp], ... % N.B.  XYZC Data must have at least 2 cols
                'YData', [abs(Jtot) abs(Jtot)], ...
                'ZData', [t' t'], ...
                'CData', [t' t'], ...
                'FaceColor', 'none', ...
                'EdgeColor', 'interp', ...
                'Marker', 'none', 'LineWidth', 1);
            s1 = gca;
            xlabel('Vapp [V]');
            ylabel('|J| [A cm^{-2}]');

            if ylogon
                set(s1, 'YScale', 'log');
            else
                set(s1, 'YScale', 'linear');
            end

            hold off
        end

        function xmesh(sol)
            figure(11)
            plot(sol.x)
            xlabel('Point')
            ylabel('Position [cm]')
        end

        % electrostatic potential as a function of position
        function Vx(varargin)
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            figure('Name', 'Vx');
            df_plot.x2d(sol, x, {V}, {'V'}, {'-'}, 'Electrostatic potential [V]', tarr, xrange, 0, 0);
        end

        % electrostatic potential as a function of position
        function Fx(varargin)
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            F = df_analysis.calcF(sol, "whole");
            figure('Name', 'Fx');
            df_plot.x2d(sol, x, {F}, {'F'}, {'-'}, 'Electric field [Vcm-1]', tarr, xrange, 0, 0);
        end

        function npx(varargin)
            % Carrier densities as a function of position
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);

            figure('Name', 'npx');
            df_plot.x2d(sol, x, {n, p}, {'n', 'p'}, {'-', '-'}, 'Carrier density [cm-3]', tarr, xrange, 0, 1)
        end

        function nspsx(varargin)
            % Carrier densities as a function of position
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            [~, ns, ps, ~, ~] = df_analysis.calcr(sol, "sub");

            figure(131);
            df_plot.x2d(sol, x, {n, p}, {'n', 'p'}, {'-', '-'}, ...
                'Carrier density [cm-3]', tarr, xrange, 0, 1)
            hold on
            df_plot.x2d(sol, par.x_sub, {ns, ps}, {'ns', 'ps'}, {'-.', '-.'}, ...
                'Carrier density [cm-3]', tarr, xrange, 0, 1)
        end

        function acx(varargin)
            % Ionic carrier densities as a function of position
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);

            Nani = repmat(dev.Nani, length(t), 1);
            Ncat = repmat(dev.Ncat, length(t), 1);

            figure('Name', 'acx');
            df_plot.x2d(sol, x, {a, c, Ncat, Nani}, {'anion', 'cation', 'static cation', 'static anion'}, {'-', '-', '--', '-.'}, ...
                'Ionic carrier density [cm-3]', tarr, xrange, 0, 0);
        end

        function gx(varargin)
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            [g1, g2, g] = df_analysis.calcg(sol);

            figure(15)
            df_plot.x2d(sol, par.x_sub, {g1, g2, g}, {'g1', 'g2', 'g total'}, ...
                {'-', '-', '-'}, 'Generation rate [cm-3s-1]', tarr, xrange, 0, 0);
        end

        function gxt(sol)
            % Carrier densities as a function of position
            par = sol.par;
            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);
            [~, ~, g] = df_analysis.calcg(sol);
            xnm = par.x_sub * 1e7;

            figure(16)
            surf(xnm, t, g)
            xlabel('Position [cm]')
            ylabel('Time [s]')
            zlabel('Generation rate [cm^{-3}s^{-1}]')
        end

        function rx(varargin)
            % Recombination rates as a function of position
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            x_sub = par.x_sub;
            r = df_analysis.calcr(sol, "sub");

            figure('Name', 'rx')
            df_plot.x2d(sol, x_sub, {r.btb, r.srh, r.vsr, r.tot}, {'rbtb', 'rsrh', 'rvsr', 'rtot'}, ...
                {'-', '-', '-', '-'}, 'Recombination rate [cm-3s-1]', tarr, xrange, 0, 0);
        end

        function rsrhx(varargin)
            % Recombination rates as a function of position
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            x = par.x_sub;
            r = df_analysis.calcr(sol, "sub");

            figure(171)
            df_plot.x2d(sol, x, {r.srh}, {''}, ...
                {'-'}, 'SRH recombination rate [cm-3s-1]', tarr, xrange, 0, 1);
        end

        function rvsrx(varargin)
            % Recombination rates as a function of position
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            x_sub = par.x_sub;
            r = df_analysis.calcr(sol, "sub");

            figure(171)
            df_plot.x2d(sol, x_sub, {r.vsr}, {''}, ...
                {'-'}, 'Surface recombination rate [cm-3s-1]', tarr, xrange, 0, 1);
        end

        function JrecVapp(JV, option)
            % Plots recombination currents for JV

            % JV - a solution from doJV
            % OPTION - 1 = dark only, 2 = light only, 3 = dark & light
            % JV is a structure containing dark and illuminated JVs

            if option == 1 || option == 3
                J.dk.f = df_analysis.calcJ(JV.dk.f);
                Vapp.dk.f = df_analysis.calcVapp(JV.dk.f);
                J.dk.r = df_analysis.calcJ(JV.dk.r);
                Vapp.dk.r = df_analysis.calcVapp(JV.dk.r);

                figure(13)
                plot(Vapp.dk.f, J.dk.f.tot(:, end), '--', Vapp.dk.r, J.dk.r.tot(:, end));
                hold on
            end

            if option == 2 || option == 3
                solf = JV.ill.f;
                solr = JV.ill.r;
                par = solf.par;
                pcum0 = par.pcum0;

                J.ill.f = df_analysis.calcJ(JV.ill.f);
                Vapp.ill.f = df_analysis.calcVapp(JV.ill.f);
                J.ill.r = df_analysis.calcJ(JV.ill.r);
                Vapp.ill.r = df_analysis.calcVapp(JV.ill.r);

                r_f = df_analysis.calcr(JV.ill.f, "whole");
                Jrec_btb_f = JV.ill.f.par.e * trapz(JV.ill.f.x, r_f.btb, 2);
                Jrec_srhint_f = JV.ill.f.par.e * trapz(JV.ill.f.x(pcum0(2) + 1:pcum0(3)), r_f.srh(:, pcum0(2) + 1:pcum0(3)), 2) ...
                    +JV.ill.f.par.e * trapz(JV.ill.f.x(pcum0(4) + 1:pcum0(5)), r_f.srh(:, pcum0(4) + 1:pcum0(5)), 2);
                Jrec_srhbulk_f = JV.ill.f.par.e * trapz(JV.ill.f.x(pcum0(3) + 1:pcum0(4)), r_f.srh(:, pcum0(3) + 1:pcum0(4)), 2);
                Jrec_tot_f = JV.ill.f.par.e * trapz(JV.ill.f.x, r_f.tot, 2);

                r_rev = df_analysis.calcr(JV.ill.r, "whole");
                Jrec_btb_r = JV.ill.f.par.e * trapz(JV.ill.r.x, r_rev.btb, 2);
                Jrec_srhint_r = JV.ill.r.par.e * trapz(JV.ill.r.x(pcum0(2) + 1:pcum0(3)), r_rev.srh(:, pcum0(2) + 1:pcum0(3)), 2) ...
                    +JV.ill.r.par.e * trapz(JV.ill.r.x(pcum0(4) + 1:pcum0(5)), r_rev.srh(:, pcum0(4) + 1:pcum0(5)), 2);
                Jrec_srhbulk_r = JV.ill.r.par.e * trapz(JV.ill.r.x(pcum0(3) + 1:pcum0(4)), r_rev.srh(:, pcum0(3) + 1:pcum0(4)), 2);
                Jrec_tot_r = JV.ill.r.par.e * trapz(JV.ill.r.x, r_rev.tot, 2);

                cc = lines(4);

                figure(13)
                plot(Vapp.ill.f, J.ill.f.tot(:, end), '--', 'Color', cc(1, :));
                hold on
                plot(Vapp.ill.r, J.ill.r.tot(:, end), 'Color', cc(1, :));
                % Recombination currents
                plot(Vapp.ill.f, Jrec_btb_f, '--', 'Color', cc(2, :));
                plot(Vapp.ill.r, Jrec_btb_r, 'Color', cc(2, :));
                plot(Vapp.ill.f, Jrec_srhint_f, '--', 'Color', cc(3, :));
                plot(Vapp.ill.r, Jrec_srhint_r, 'Color', cc(3, :));
                plot(Vapp.ill.f, Jrec_srhbulk_f, '--', 'Color', cc(4, :));
                plot(Vapp.ill.r, Jrec_srhbulk_r, 'Color', cc(4, :));
            end

            figure(13)
            ylim([-30e-3, 10e-3]);
            xlabel('Applied voltage [V]')
            ylabel('Current density [Acm-2]');
            hold off
            legend('Illumated for', 'Illumated rev', 'Jrec,btb for', 'Jrec, btb rev' ...
                , 'Jrec,srh-int for', 'Jrec,srh-int rev', 'Jrec,srh-bulk for', 'Jrec,srh-bulk rev')
        end

        function Ft(sol, xpos)
            % Absolute field strength F as a function of time at point
            % position XPOS
            [~, t, xmesh, ~, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);
            ppos = getpointpos(xpos, xmesh);

            F = df_analysis.calcF(sol, "whole");

            figure(14)
            plot(t, F(:, ppos))
            xlabel('Time [s]')
            ylabel(['Electric Field at pos x = ', num2str(round(xpos * 1e7)), 'nm [Vcm-1]'])
        end

        function sigmat(sol)
            % Plot the integrated space charge density [cm-2] as a function of time
            sigma = df_analysis.calcsigma(sol);
            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);
            figure(15)
            plot(t, sigma)
            xlabel('Time [s]')
            ylabel('sigma [C cm-2]')
        end

        function Qt(sol, x1, x2)
            % Plot the integrated space charge density in Coulombs [Ccm-2] as a function of time
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);

            p1 = find(x <= x1);
            p1 = p1(end);
            p2 = find(x <= x2);
            p2 = p2(end);

            rho = df_analysis.calcrho(sol, "whole");
            Q = par.e * trapz(x(p1:p2), rho(:, p1:p2), 2); % net charge

            figure('Name', 'Q v.s. time');
            plot(t, Q)
            xlabel('Time [s]')
            ylabel('Charge [C cm-2]')
            xlim([t(1), t(end)])
        end

        function Qionict(sol, x1, x2)
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            p1 = find(x <= x1); p1 = p1(end);
            p2 = find(x <= x2); p2 = p2(end);
            rho_ionic = df_analysis.rhoIonicCalc(sol, "whole");
            Q_ionic = par.e * trapz(x(p1:p2), rho_ionic(:, p1:p2), 2); % net charge

            figure('Name', 'Q(ionic) v.s. time');
            yyaxis left
            plot(t, Q_ionic)
            xlabel('Time [s]')
            ylabel('Charge [C cm-2]')
            xlim([t(1), t(end)])
        end

        % integrated charge density as a function of applied voltage
        function QVapp(sol, x1, x2)
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            p1 = find(x <= x1);
            p1 = p1(end);
            p2 = find(x <= x2);
            p2 = p2(end);

            rho = df_analysis.calcrho(sol, "whole");
            Vapp = df_analysis.calcVapp(sol);
            Q = par.e * trapz(x(p1:p2), rho(:, p1:p2), 2);

            figure(17)
            plot(Vapp, Q)
            xlabel('Vapp [V]')
            ylabel('Charge [C cm-2]')

            if Vapp(1) ~= Vapp(end)
                xlim([Vapp(1), Vapp(end)])
            end

        end

        function rhox(varargin)
            % Volumetric charge density (rho) as a funciton of position
            % A time array can be used as a second input argument
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            rho = df_analysis.calcrho(sol, "whole");

            figure(19)
            df_plot.x2d(sol, x, {rho}, {'\rho'}, {'-'}, 'Charge density [cm-3]', tarr, xrange, 0, 0);
        end

        function deltarhox(varargin)
            % the change in volumetric charge density (rho) as a funciton of position
            % a time array can be used as a second input argument
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            rho = df_analysis.calcrho(sol, "whole");
            deltarho = rho - rho(1, :);

            figure(20)
            df_plot.x2d(sol, x, {deltarho}, {'\Delta \rho'}, {'-'}, 'Delta charge density [cm-3]', tarr, xrange, 0, 0);
        end

        function rhoxFxVx(varargin)
            % Three panel figure:
            % Volumetric charge density (rho), Field and potential as a funciton of position
            % A time array can be used as a second input argument
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);

            rho = df_analysis.calcrho(sol, "whole");
            F = df_analysis.calcF(sol, "whole");

            figure(21)
            subplot(3, 1, 1)
            df_plot.x2d(sol, x, {rho}, {'\rho'}, {'-'}, 'Charge density [cm-3]', tarr, xrange, 0, 0);
            subplot(3, 1, 2)
            df_plot.x2d(sol, x, {F}, {'F'}, {'-'}, 'Electric field [Vcm-1]', tarr, xrange, 0, 0);
            subplot(3, 1, 3)
            df_plot.x2d(sol, x, {V}, {'V'}, {'-'}, 'Electrostatic potential [V]', tarr, xrange, 0, 0);
        end

        function rhoxVx(varargin)
            % Three panel figure:
            % Volumetric charge density (rho), Field and potential as a funciton of position
            % A time array can be used as a second input argument
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);

            rho = df_analysis.calcrho(sol, "whole");

            figure(211)
            subplot(2, 1, 1)
            df_plot.x2d(sol, x, {rho}, {'\rho'}, {'-'}, 'Charge density [cm-3]', tarr, xrange, 0, 0);
            subplot(2, 1, 2)
            df_plot.x2d(sol, x, {-V}, {'V'}, {'-'}, '-Electrostatic potential [V]', tarr, xrange, 0, 0);
        end

        % energy Level diagram, and charge densities plotter
        function ELx(varargin)
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = df_analysis.calcEnergies(sol);

            figure('Name', 'ELx');
            df_plot.x2d(sol, x, {Efn, Efp, Ecb, Evb}, {'E_{fn}', 'E_{fp}', 'E_{CB}', 'E_{VB}'}, ...
                {'--', '--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0)
        end

        % energy Level diagram, and charge densities plotter
        function ELx_uncontacted(varargin)
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            Ecb = repmat(dev.Phi_EA, length(t), 1);
            Evb = repmat(dev.Phi_IP, length(t), 1);
            EF0 = repmat(dev.EF0_zerointerface, length(t), 1);

            figure(22);
            df_plot.x2d(sol, x, {EF0, Ecb, Evb}, {'E_{F0}', 'E_{CB}', 'E_{VB}'}, ...
                {'--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0)
        end

        % energy Level diagram, and charge densities plotter
        function ELnpx(varargin)
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = df_analysis.calcEnergies(sol);

            figure('Name', 'ELnpx')
            subplot(2, 1, 1);
            df_plot.x2d(sol, x, {Efn, Efp, Ecb, Evb}, {'E_{fn}', 'E_{fp}', 'E_{CB}', 'E_{VB}'}, ...
                {'--', '--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0);
            subplot(2, 1, 2);
            df_plot.x2d(sol, x, {n, p}, {'n', 'p'}, {'-', '-'}, 'El carrier density [cm-3]', tarr, xrange, 0, 1);
        end

        % energy Level diagram, and charge densities plotter
        function ELxnpxacx(varargin)
            % SOL = the solution structure
            % TARR = An array containing the times that you wish to plot
            % XRANGE = 2 element array with [xmin, xmax]
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = df_analysis.calcEnergies(sol);
            NA = repmat(dev.NA, length(t), 1);
            ND = repmat(dev.ND, length(t), 1);

            Nani = repmat(dev.Nani, length(t), 1);
            Ncat = repmat(dev.Ncat, length(t), 1);

            figure('Name', 'ELxnpxacx');
            subplot(3, 1, 1);
            df_plot.x2d(sol, x, {Efn, Efp, Ecb, Evb}, {'E_{fn}', 'E_{fp}', 'E_{CB}', 'E_{VB}'}, {'--', '--', '-', '-'}, 'Energy [eV]', tarr, xrange, 0, 0);
            subplot(3, 1, 2);
            df_plot.x2d(sol, x, {n, p}, {'electrons, \it{n}', 'holes, \it{p}'}, {'-', '-'}, 'Density [cm-3]', tarr, xrange, 0, 1);
            subplot(3, 1, 3);
            df_plot.x2d(sol, x, {a, c, Ncat, Nani}, {'anion', 'cation', 'static cation', 'static anion'}, {'-', '-', '--', '-.'}, ...
                'Density [cm-3]', tarr, xrange, 0, 0);
        end

        % potential and ionic charges as a function of position
        function Vxacx(varargin)
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            Nani = repmat(dev.Nani, length(t), 1);
            Ncat = repmat(dev.Ncat, length(t), 1);

            figure('Name', 'Vxacx')
            subplot(2, 1, 1);
            df_plot.x2d(sol, x, {V}, {'V'}, ...
                {'-'}, 'Electro. potential [V]', tarr, xrange, 0, 0);
            subplot(2, 1, 2);
            df_plot.x2d(sol, x, {a, c, Ncat, Nani}, {'anion', 'cation', 'static cation', 'static anion'}, {'-', '-', '--', '-.'}, ...
                'Ionic carrier density [cm-3]', tarr, xrange, 0, 0);
        end

        % electrostatic potential as a function of position
        function Vionxacx(varargin)
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            Vion = df_analysis.calcVion(sol);
            Vel = V - Vion;

            figure('Name', 'Vionxacx')
            subplot(2, 1, 1);
            df_plot.x2d(sol, x, {V, Vion, Vel}, {'V', 'Vion', 'Vel'}, ...
                {'--', '.', '-'}, 'Electro. potential [V]', tarr, xrange, 0, 0);
            subplot(2, 1, 2);
            df_plot.x2d(sol, x, {a, c}, {'a', 'c'}, ...
                {'-', '-'}, 'Ionic carrier density [cm-3]', tarr, xrange, 0, 0);
        end

        function Fiont(sol, xpos)
            % field contribution from ionic charge FION as a function of time at position XPOS
            [~, t, xmesh, ~, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);
            ppos = getpointpos(xpos, xmesh);
            Fion = df_analysis.calcFion(sol);

            figure(26)
            plot(t, Fion(:, ppos))
            xlabel('Time')
            ylabel('Ion field [Vcm-1]')
        end

        function rec_zone(sol)
            dev = sol.par.dev_sub;
            x = sol.par.x_sub;

            figure(27)
            plot(x, dev.int_switch, x, dev.srh_zone, x, dev.vsr_zone)
            xlabel('Position [nm]')
            ylabel('norm')
            legend('interface', 'SRH zone', 'VSR zone')
        end

        function alpha0beta0(sol)
            dev = sol.par.dev_sub;
            x = sol.par.x_sub;

            figure(28)
            plot(x, dev.alpha0, x, dev.beta0, x, dev.alpha0_xn, x, dev.beta0_xp)
            xlabel('Position [nm]')
            ylabel('')
            legend('alpha0', 'beta0', 'alpha0-xn', 'beta0-xp')

        end

        function colourblocks(sol, yrange)
            par = sol.par;
            dcum0 = par.dcum0 * 1e7; % Convert to nm

            for i = 1:length(dcum0) - 1
                v = [dcum0(i) yrange(2); dcum0(i + 1) yrange(2); dcum0(i + 1) yrange(1); dcum0(i) yrange(1)]; % vertices position
                f = [1 2 3 4]; % Faces

                if length(par.layer_colour) == length(dcum0) - 1
                    j = i;
                else
                    j = i - ((length(par.layer_colour) - 1) * floor(i / length(par.layer_colour)));
                end

                colour = par.layer_colour(j, :);
                patch('Faces', f, 'Vertices', v, 'FaceColor', colour, 'EdgeColor', 'none'); %,'HandleVisibility', 'off')
            end

            hold on
        end

        function [sol, tarr, pointtype, xrange] = sortarg(args)

            if length(args) == 1
                sol = args{1};
                tarr = sol.t(size(sol.u, 1));
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)] * 1e7; % converts to nm
            elseif length(args) == 2
                sol = args{1};
                tarr = args{2};
                pointtype = 't';
                xrange = [sol.x(1), sol.x(end)] * 1e7; % converts to nm
            elseif length(args) == 3
                sol = args{1};
                tarr = args{2};
                xrange = args{3};
                pointtype = 't';
            end

        end

        function x2d(sol, xmesh, variables, legstr, linestyle, ylab, tarr, xrange, logx, logy)
            % SOL = solution structure
            % VARIABLES is an array containing the variables for plotting
            % LEGSTR is the legend string
            % YLAB = y-axis label
            % TARR- array of times
            % XRANGE - limits of the plot as a two element vector
            % LOGX, LOGY - switches for log axes
            ax = gca;

            if ishold(ax) == 0
                cla(ax); % Clear current axis if held
            end

            par = sol.par;
            xnm = xmesh * 1e7;

            vmin = min(min(cell2mat(variables)));
            vmax = max(max(cell2mat(variables)));

            if vmin == 0 && vmax == 0
                vmin = -1;
                vmax = 1;
            end

            vrange = vmax - vmin;

            if isempty(findobj(ax, 'Type', 'patch'))

                switch logy
                    case 0
                        df_plot.colourblocks(sol, [vmin - (vrange * 0.2), vmax + (vrange * 0.2)]);
                    case 1
                        df_plot.colourblocks(sol, [0.1 * vmin, 10 * vmax]);
                end

            end

            vmin_tarr = zeros(length(tarr), length(variables));
            vmax_tarr = zeros(length(tarr), length(variables));
            h = zeros(1, length(variables));

            for i = 1:length(tarr)
                % find the time
                p1 = find(sol.t <= tarr(i));
                p1 = p1(end);

                for jj = 1:length(variables)
                    vtemp = variables{jj};

                    vmin_tarr(i, jj) = min(vtemp(p1, :));
                    vmax_tarr(i, jj) = max(vtemp(p1, :));

                    h(i, jj) = plot(xnm, variables{jj}(p1, :), char(linestyle(jj)));
                    hold on
                end

            end

            xlabel('Position [nm]')
            ylabel(ylab)

            if logy == 1
                set(gca, 'YScale', 'log');
            end

            if logx == 1
                set(gca, 'XScale', 'log');
            end

            if length(variables) == 1
                mystr = [];

                for i = 1:length(tarr)
                    mystr = [mystr, string(['t = ', num2str(tarr(i)), ' s'])];
                end

                lgd = legend(h, mystr);
            else
                lgd = legend(h(1, :), legstr);
            end

            lgd.FontSize = 12;
            xlim([xrange(1), xrange(2)])
            ymin = min(min(vmin_tarr));
            ymax = max(max(vmax_tarr));
            yrange = ymax - ymin;

            if ymin == 0 && ymax == 0
            else

                switch logy
                    case 0

                        if yrange == 0
                            ylim([ymin * 0.9, ymax * 1.1]);
                        else
                            ylim([ymin - (yrange * 0.2), ymax + (yrange * 0.2)]);
                        end

                    case 1
                        ylim([0.1 * ymin, 10 * ymax])
                end

            end

            set(gca, 'Layer', 'top')
            box on
            hold off
        end

        % ADD-ONs: cycle-to-cycle (C2C) characterisation
        function c2c(sol, xpos, area_coeff, cycle, use_abs)

            J = df_analysis.calcJ(sol); % current flux (can select different species)
            Vapp = df_analysis.calcVapp(sol); % applied voltage
            t = sol.t;
            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);
            data_pnts_per_cyc = length(Vapp) / (cycle);

            legend_handle = []; % legend handles
            legend_label = cell(1, cycle);

            figure('Name', 'Cycle-to-cycle Variability');
            hold on;

            for i = 1:cycle

                % round (), make sure idx is an integer
                % (1 : data_pnts_per_cyc), current data points array
                % e.g., round((1 - 1) * (data_pnts_per_cyc), 1st Cycle starts from 0
                %       round((2 - 1) * (data_pnts_per_cyc), 2nd Cycle starts after 1st Cycle

                cycle_idx = round(round(i - 1) * (data_pnts_per_cyc) + (1:data_pnts_per_cyc));

                if use_abs == 1
                    plots = plot(Vapp(cycle_idx), abs(J.tot(cycle_idx, ppos)) * area_coeff, 'LineWidth', 0.75);
                else
                    plots = plot(Vapp(cycle_idx), (J.tot(cycle_idx, ppos) * area_coeff), 'LineWidth', 0.75);
                end

                legend_handle = [legend_handle, plots];
                legend_label{i} = ['C' num2str(i)];
            end

            hold off;
            box on;
            xlabel('Applied Voltage, Vapp [V]');
            ylabel('Current (A)');
            legend(legend_handle, legend_label, 'Location', 'northwest', 'FontSize', 10);
            % yscale log;
        end

        % ADD-ONs: Current-Vapp track
        function CVappTrack(sol, xpos, current_type, area_coeff, absolutely, logarithmically)

            % data processing
            J = df_analysis.calcJ(sol);
            Vapp = df_analysis.calcVapp(sol);
            % mesh
            t = sol.t; % time mesh
            xmesh = sol.x; % spatial mesh
            ppos = getpointpos(xpos, xmesh);
            % plotting
            figure('Name', 'Current-Voltage Plot')

            % select the current type plot
            switch current_type

                case "all"

                    if absolutely == 1
                        plot(Vapp, abs(J.n(:, ppos)) * area_coeff, ...
                            Vapp, abs(J.p(:, ppos)) * area_coeff, ...
                            Vapp, abs(J.c(:, ppos)) * area_coeff, ...
                            Vapp, abs(J.a(:, ppos)) * area_coeff, ...
                            Vapp, abs(J.disp(:, ppos)) * area_coeff, ...
                            Vapp, abs(J.tot(:, ppos)) * area_coeff ...
                        );
                    else
                        plot(Vapp, J.n(:, ppos) * area_coeff, ...
                            Vapp, J.p(:, ppos) * area_coeff, ...
                            Vapp, J.c(:, ppos) * area_coeff, ...
                            Vapp, J.a(:, ppos) * area_coeff, ...
                            Vapp, J.disp(:, ppos) * area_coeff, ...
                            Vapp, J.tot(:, ppos) * area_coeff ...
                        );
                    end

                    legend('In', 'Ip', 'Ic', 'Ia', 'Idisp', 'Itotal');

                case "total"
                    % * for total current plot
                    vpnts = length(Vapp);
                    ref = floor(vpnts / 4);

                    part_first = 1:ref;
                    part_second = (ref + 1):(2 * ref);
                    part_third = (2 * ref + 1):(3 * ref);
                    part_fourth = (3 * ref + 1):vpnts;

                    J_total = J.tot(:, ppos);

                    % the scan track only works for the symmetric voltage applied
                    % e.g., (1, -1) or (-1, 1)
                    % not work for asymmetric voltage applied
                    % e.g., (1.2, -0.8)

                    hold on;

                    if absolutely == 1
                        plot(Vapp(part_first), abs(J_total(part_first)) * area_coeff, 'r.');
                        plot(Vapp(part_second), abs(J_total(part_second)) * area_coeff, 'y.');
                        plot(Vapp(part_third), abs(J_total(part_third)) * area_coeff, 'g.');
                        plot(Vapp(part_fourth), abs(J_total(part_fourth)) * area_coeff, 'b.');
                    else
                        plot(Vapp(part_first), J_total(part_first) * area_coeff, 'r.');
                        plot(Vapp(part_second), J_total(part_second) * area_coeff, 'y.');
                        plot(Vapp(part_third), J_total(part_third) * area_coeff, 'g.');
                        plot(Vapp(part_fourth), J_total(part_fourth) * area_coeff, 'b.');
                    end

                    hold off;
                    legend('1st Part', '2nd Part', '3rd Part', '4th Part');

                case "electron"
                    plot(Vapp, J.n(:, ppos) * area_coeff);
                case "hole"
                    plot(Vapp, J.p(:, ppos) * area_coeff);
                case "anion"
                    plot(Vapp, J.a(:, ppos) * area_coeff);
                case "cation"
                    plot(Vapp, J.c(:, ppos) * area_coeff);
                case "disp"
                    plot(Vapp, J.disp(:, ppos) * area_coeff);
                otherwise

                    if absolutely == 1
                        plot(Vapp, abs(J.tot(:, ppos)) * area_coeff);
                    else
                        plot(Vapp, J.tot(:, ppos) * area_coeff);
                    end

            end

            xlabel('Applied Voltage, Vapp [V]');
            ylabel('Current Density, J [A cm^{-2}]');

            if logarithmically == 1 % plot y-axis in logarithm or not
                yscale log;
            end

        end

        % ADD-ONs: device-to-device (D2D) characterisation with one variable
        function d2d_var(sols, xpos)

            figure('Name', 'd2d');
            hold on;

            for i = 1:length(sols)

                J = df_analysis.calcJ(sols{i});
                Vapp = df_analysis.calcVapp(sols{i});
                t = sols{i}.t;
                xmesh = sols{i}.x;
                ppos = getpointpos(xpos, xmesh);

                plot(Vapp, abs(J.tot(:, ppos)), 'DisplayName', ['D' num2str(i)], 'LineWidth', 0.5);
            end

            hold off;
            xlabel('Applied Voltage, Vapp [V]');
            ylabel('Current Density, J [A cm^{-2}]');
            legend('Location', 'northwest', 'FontSize', 10);
            yscale log; % logrithm setting
        end

        % ADD-ONs: D2D characterisation with two variables
        function d2d_vars(sols, xpos, var_1, var_2)

            figure('Name', ['D2D (two variables)']);
            hold on;

            for i = 1:size(sols, 1)

                for j = 1:size(sols, 2)
                    xmesh = sols{i, j}.x;
                    ppos = getpointpos(xpos, xmesh);
                    J = df_analysis.calcJ(sols{i, j});
                    Vapp = df_analysis.calcVapp(sols{i, j});
                    t = sols{i, j}.t;

                    plot(Vapp, J.tot(:, ppos), ...
                        'LineWidth', 0.5, ...
                        'DisplayName', ['var[1] = ', num2str(var_1(i)), ', var[2] = ', num2str(var_2(j))]);

                end

            end

            xlabel('Applied Voltage (V)');
            ylabel('Current Density (A/cm^2)');
            legend('Location', 'northwest', 'FontSize', 10);
            hold off;
        end

        % ADD-ONs: current (species) v.s. time
        function Jt_species(sol, xpos)
            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);
            [J, j, xmesh] = df_analysis.calcJ(sol);
            ppos = getpointpos(xpos, xmesh);

            figure('Name', 'J-t (species)');
            plot(t, J.n(:, ppos), ...
                t, J.p(:, ppos), ...
                t, J.a(:, ppos), ...
                t, J.c(:, ppos), ...
                t, J.disp(:, ppos), ...
                'LineWidth', 0.15);
            legend('Jn', 'Jp', 'Ja', 'Jc', 'Jdisp');
            xlabel('time [s]');
            ylabel('J [A cm^{-2}]');
            set(legend, 'FontSize', 16);
            set(legend, 'EdgeColor', [1 1 1]);
        end

        % ADD-ONs: overpoential v.s. species current flux
        function etaJ(sol, calc)
            x_sol = sol.x;
            t_sol = sol.t;
            t_calc = calc.t; % t from solution & boundary condition calculation
            eta_calc = calc.eta; % overpotential
            [t_sparsif, sparsif_idx] = unique(t_calc, 'stable'); % sparse the data (keep the first occurrence)
            eta_sparsif = eta_calc(sparsif_idx);
            eta_sol = interp1(t_sparsif, eta_sparsif, t_sol);
            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = df_analysis.splitsol(sol);
            [J, j, xmesh] = df_analysis.calcJ(sol);

            figure('Name', 'eta - current');
            plot(eta_sol, J.tot, LineWidth = 0.5); % J.tot; J.a; J.c; J_p; J_n
            xlabel('Overpotential, eta');
            ylabel('Current Flux');
        end

        % ADD-ONs: charge carrier sepcies
        function npxacx(varargin)
            [sol, tarr, pointtype, xrange] = df_plot.sortarg(varargin);
            [u, t, x, par, dev, n, p, a, c, V] = df_analysis.splitsol(sol);
            [Ecb, Evb, Efn, Efp] = df_analysis.calcEnergies(sol);

            NA = repmat(dev.NA, length(t), 1);
            ND = repmat(dev.ND, length(t), 1);
            Nani = repmat(dev.Nani, length(t), 1);
            Ncat = repmat(dev.Ncat, length(t), 1);

            figure('Name', 'ELxnpxacx');
            subplot(2, 1, 1);
            df_plot.x2d(sol, x, {n, p}, {'electrons, \it{n}', 'holes, \it{p}'}, {'-', '-'}, 'Density [cm-3]', tarr, xrange, 0, 1);
            subplot(2, 1, 2);
            df_plot.x2d(sol, x, {a, c, Ncat, Nani}, {'anion', 'cation', 'static cation', 'static anion'}, {'-', '-', '--', '-.'}, ...
                'Density [cm-3]', tarr, xrange, 0, 0);
        end

    end

end
