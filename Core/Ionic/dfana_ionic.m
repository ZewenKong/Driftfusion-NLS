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

            j_bv = j0 * (exp((1 - alpha_e) * F * (E - E_eq) / (R * T)) - exp((- alpha_e) * F * (E - E_eq) / (R * T)));
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

        %% - - - - - - - - - - NORMALISED MOBILITY INVESTIGATION - - - - - - - - - -

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
