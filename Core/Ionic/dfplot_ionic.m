classdef dfplot_ionic

    % Driftfusion-mem dfplot content:
    % dfplot.JXT

    %% - - - - - - - - - - CODE START - - - - - - - - - -

    methods (Static)

        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        function c2c(J, Vapp, ppos, cycle)

            % J, current flux (can select different species)
            % Vapp, applied voltage
            % ppos, position point (xpos)
            % cycle

            figure('Name', 'Cycle-to-cycle Variability');
            hold on;

            data_pnts_per_cyc = length(Vapp) / (cycle);

            % - - - - - - - - - - handles
            legend_handle = [];
            legend_label = cell(1, cycle);

            for i = 1:cycle

                % round (), make sure idx is an integer
                % (1 : data_pnts_per_cyc), current data points array
                % e.g., round((1 - 1) * (data_pnts_per_cyc), 1st Cycle starts from 0
                %       round((2 - 1) * (data_pnts_per_cyc), 2nd Cycle starts after 1st Cycle

                cycle_idx = round(round(i - 1) * (data_pnts_per_cyc) + (1:data_pnts_per_cyc));
                plots = plot(Vapp(cycle_idx), J.tot(cycle_idx, ppos), 'LineWidth', 0.75);
                legend_handle = [legend_handle, plots];
                legend_label{i} = ['C' num2str(i)];
            end

            hold off;
            box on;

            xlabel('Applied Voltage, Vapp [V]');
            ylabel('Current Density, J [A cm^{-2}]');
            legend(legend_handle, legend_label, 'Location', 'northwest', 'FontSize', 10);
        end

        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        function d2d(sols, var_1, var_2)

            var_1_name = inputname(2);
            var_2_name = inputname(3);

            figure('Name', ['D2D (', var_1_name, ' & ', var_2_name, ')']);
            hold on;

            for i = 1:size(sols, 1)

                for j = 1:size(sols, 2)
                    xmesh = sols{i, j}.x;
                    xpos = 0;
                    ppos = getpointpos(xpos, xmesh);

                    J = dfana.calcJ(sols{i, j});
                    Vapp = dfana.calcVapp(sols{i, j});
                    t = sols{i, j}.t;

                    plot(Vapp, J.tot(:, end), ...
                        'LineWidth', 0.5, ...
                        'DisplayName', ['j0 = ', num2str(var_1(i)), ' B\_ionic = ', num2str(var_2(j))]);

                end

            end

            xlabel('Applied Voltage (V)');
            ylabel('Current Density (A/cm^2)');
            legend('Location', 'northwest', 'FontSize', 10);
            hold off;
        end

        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        function Jt_species(sol, xpos)

            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = dfana.splitsol(sol);

            [J, j, xmesh] = dfana.calcJ(sol);
            ppos = getpointpos(xpos, xmesh);

            figure(101);
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

        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        function eta_flux(sol, calc)
            %
            % overpotential v.s. species current flux
            %
            x_sol = sol.x;
            t_sol = sol.t; t_calc = calc.t; % t from solution & boundary condition calculation
            eta_calc = calc.eta; % overpotential

            [t_sparsif, sparsif_idx] = unique(t_calc, 'stable'); % sparse the data (keep the first occurrence)
            eta_sparsif = eta_calc(sparsif_idx);
            eta_sol = interp1(t_sparsif, eta_sparsif, t_sol);
            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = dfana.splitsol(sol);
            [J, j, xmesh] = dfana.calcJ(sol);

            figure('Name', 'Cycle-to-cycle Variability');

            plot(eta_sol, J.tot, LineWidth = 0.5); % J.tot; J.a; J.c; J_p; J_n

            xlabel('Overpotential, eta');
            ylabel('Current Flux');

            % - - - - - - - - - - Adds-on (time, position)
            %
            % - - - - - - - - - - v.s. time
            % plot(t_calc, eta_calc, '.');

            % - - - - - - - - - - v.s. position
            % x_sparsif = linspace(x_sol(1), x_sol(end), length(eta_sol));
            % plot(x_sparsif, eta_sol);
        end

        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        function Qt_ionic(sol, x1, x2)

            [u, t, x, par, dev, n, p, a, c, V] = dfana.splitsol(sol);

            p1 = find(x <= x1); p1 = p1(end);
            p2 = find(x <= x2); p2 = p2(end);

            rho_ionic = dfana_ionic.calcrho_ionic(sol, "whole");
            Q_ionic = par.e * trapz(x(p1:p2), rho_ionic(:, p1:p2), 2); % net charge

            figure('Name', 'Q(ionic) v.s. time');
            yyaxis left
            plot(t, Q_ionic)
            xlabel('Time [s]')
            ylabel('Charge [C cm-2]')
            xlim([t(1), t(end)])

            % add second y-axis for Vapp vs time
            yyaxis right
            Vapp = dfana.calcVapp(sol);
            plot(t, Vapp);
        end

        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    end

end
