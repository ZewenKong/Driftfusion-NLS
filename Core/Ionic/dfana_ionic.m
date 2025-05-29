classdef dfana_ionic

    methods (Static)

        function a = denstochemact(density) % transfers charge carrier density to chemical activity
            A = 6.022e23; % avogadro const.
            a = (density * 10e3) / A; % (cm-3 to dm-3 (SI unit)) / avogadro const.

            % from wikipedia: for a given dissolved species, its chemical activity (a)
            % is the product of its activity coefficient (γ) by its molar (mol/L solution),
            % or molal (mol/kg water), concentration (C): a = γ C.
        end

        function E_eq = nernst(a, E_st, R, T, F, z) % Nernst equation calculation

            % a, the chemical activity
            % E_st, the standard potential (Ag/AgI = 0.152)
            % R, Universal ideal gas constant
            % T, Temperature
            % F, Faraday constant
            % z, number of charge transfer (1)

            E_eq = E_st - (((R * T) / F) / z) * log(a);
        end

        function j_bv = butlervolmer(j0, alpha_e, R, T, F, E, E_eq) % Butler-Volmer calculation (ionic current flux)

            % j0, exchange current density (no bias)
            % alpha_e, electrode (anodic/cathodic) charge transfer coefficient, normally 1/2
            % E is the electrode potential on the electrochemical scale
            % E_eq, is the equilibrium potential calculated by Nernst equation on the electrochemical scale

            j_bv =- j0 * (exp((1 - alpha_e) * F * (E - E_eq) / (R * T)) - exp((- alpha_e) * F * (E - E_eq) / (R * T))); % negative sign accounts for current direction being opposite to driftfusion direction sign

        end

        function save_calcdata(varargin)
            persistent calc;

            if nargin == 1 && ischar(varargin{1}) && strcmp(varargin{1}, 'reset')
                calc = []; % reset the persistent variable
                disp('dfana_ionic.m - func save_calcdata: persistent calc has been reset');
                disp('-');
                return;
            end

            if isempty(calc)
                calc = struct();
            end

            for i = 1:2:nargin
                field_name = varargin{i}; % field name
                field_value = varargin{i + 1}; % field value

                if ~isfield(calc, field_name)
                    calc.(field_name) = []; % initialise the field if it doesn't exist
                end

                calc.(field_name)(end + 1, 1) = field_value; % append the value to the field
            end

            assignin('base', 'calc', calc);
        end

        function rho_ionic = calcrho_ionic(sol, mesh_option)

            [u, t, x, par, dev_in, n_whole, p_whole, a_whole, c_whole, V_whole] = dfana.splitsol(sol);

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

        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        % after C2C simulation, this file is used to investigate the electronic current flux

        function [J_max, J_corr, idx_J_max, idx_J_corr] = JmaxANDcorr(sol, xpos)

            % function to find max J and its corresponding value

            J = dfana.calcJ(sol);
            Vapp = dfana.calcVapp(sol);

            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);

            J_temp = J.tot(:, ppos);
            Vapp_temp = Vapp;

            % get max J and its idx
            [J_max, idx_J_max] = max(J_temp);

            % get max J corresponding Vapp
            V_corr_J_max = Vapp_temp(idx_J_max);

            % filter the Vapp_temp
            % find the V values smaller than / equal to V_corr_J_max
            idx_filtered_arr = find(Vapp_temp <= V_corr_J_max);

            % transfer to the max diff array
            % get the max diff value
            [max_diff, idx_max_diff] = max(diff(idx_filtered_arr));

            % plus 1 to get the index of minuend
            % which is the index of corresponding J value
            idx_J_corr = idx_filtered_arr(idx_max_diff + 1);
            J_corr = J_temp(idx_J_corr);
        end

        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        function [t_J_max, t_J_corr] = tPntOfSpJs(sol, xpos)

            % time points of special current fluxes

            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = dfana.splitsol(sol);
            J = dfana.calcJ(sol);
            Vapp = dfana.calcVapp(sol);

            xmesh = sol.x;
            ppos = getpointpos(xpos, xmesh);

            J_temp = J.tot(:, ppos);
            Vapp_temp = Vapp;

            [J_max, J_corr] = dfana_ionic.JmaxANDcorr(sol, xpos);

            idx_J_max_1 = find(J_temp == J_max);
            idx_J_corr_1 = find(J_temp == J_corr);

            t_J_max = t(idx_J_max_1);
            t_J_corr = t(idx_J_corr_1);

        end

        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        % different points investigation

        function [selsoleq_1, selsoleq_2] = spPntTosoleq(sol, xpos) % selected (max J and its corr) point to equilibrium solution

            % - - - - - - - - - - solution processing

            [~, t, ~, ~, ~, ~, ~, ~, ~, ~] = dfana.splitsol(sol);
            J = dfana.calcJ(sol);
            Vapp = dfana.calcVapp(sol);

            % - - - - - - - - - - limitation

            [J_max, J_corr, idx_J_max, idx_J_corr] = dfana_ionic.JmaxANDcorr(sol, xpos);
            [t_J_max, t_J_corr] = dfana_ionic.tPntOfSpJs(sol, xpos);
            idx_t_J_max = find(t == t_J_max);
            idx_t_J_corr = find(t == t_J_corr);

            % - - - - - - - - - - point of max J as new equilibrium solution

            selsoleq_1 = sol;
            selsoleq_1.u = sol.u(idx_t_J_max, :, :);
            selsoleq_1.t = sol.t(idx_t_J_max);
            V_corr_J_max = Vapp(idx_J_max);

            selsoleq_1.par.mobseti = 0; % * set mobseti = 0;

            [~, sol_dwell_1] = ramped_step(selsoleq_1, V_corr_J_max, 0.1, 0.1); % ramp the solution

            sol_dwell_1.par.mobseti = 1; % set mobseti back

            selsoleq_1 = sol_dwell_1;

            % - - - - - - - - - - corrsponding point of max J as new equilibrium solution

            selsoleq_2 = sol;
            selsoleq_2.u = sol.u(idx_t_J_corr, :, :);
            selsoleq_2.t = sol.t(idx_t_J_corr);
            V_corr_J_corr = Vapp(idx_J_corr);

            selsoleq_2.par.mobseti = 0;

            [~, sol_dwell_2] = ramped_step(selsoleq_2, V_corr_J_corr, 0.1, 0.1);

            sol_dwell_2.par.mobseti = 1;

            selsoleq_2 = sol_dwell_2;

        end

        % ############################################################
        % ### NORMALISED MOBILITY INVESTIGATION  #####################
        % ############################################################

        function [F_idx, maxF_idx] = maxFatIF(sol)

            [u, t, x_whole, par, dev, n, p, a, c, V] = dfana.splitsol(sol);

            F = dfana.calcF(sol, "whole"); % the potential through the device
            tot_points = cumsum(par.layer_points); % accumulate the position points
            boundary_idx = tot_points(length(par.layer_points) - 2); % points at the boundary (interface) index
            F_idx = F(:, boundary_idx); % all potential (over time) at that boundary point
            maxF_idx = max(F_idx); % max potential at that point
        end

        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

            % - - - - - - - - - - Probability term

            Ea = 0.45 * e; % activation energy (eV), approximate value
            v0 = 5e12; % phonon frequency (General value)
            c = 1; % for a 1D transport, approximately equal to unity

            % - - - - - - - - - - Equation

            P = c * v0 * exp(- (Ea / (kB * T))); % transition probability

            vD = d * P * (exp((z * e * d * F) / (2 * kB * T)) ...
                - exp(- (z * e * d * F) / (2 * kB * T)));
            vD0 = d * P * ((z * e * d * F) / (kB * T));

            mu = vD ./ F; % mobility (drfit velocity divide field potential)
            mu0 = vD0 ./ F; % mobility at low electric field (dF/kBT << 1)
            mun = mu ./ mu0;

            % The non-linear transport
            % non-linear --> longer retention time
            % just try different range, verigy w/ experimental staff
            % Plot: mun, Ea, thickness
        end

        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    end

end
