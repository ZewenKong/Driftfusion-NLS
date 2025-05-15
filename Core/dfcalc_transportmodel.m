classdef dfcalc_transportmodel

    methods (Static)
        
        function [F_idx, maxF_idx] = maxFatIF(sol)

            % find max potential at the interface
            
            [u, t, x_whole, par, dev, n, p, a, c, V] = dfana.splitsol(sol);
            F = dfana.calcF(sol, "whole");
            tot_points = cumsum(par.layer_points);
            boundary_idx = tot_points(length(par.layer_points) - 2);

            F_idx = F(:, boundary_idx);
            maxF_idx = max(F_idx);
        end

        function maxF = maxFwitinActive(sol)

            [u, t, x_whole, par, dev, n, p, a, c, V] = dfana.splitsol(sol);
            F = dfana.calcF(sol, "whole");

            tot_points = cumsum(par.layer_points);
            active_start_idx = tot_points(2) + 1;
            active_end_idx = tot_points(3);
            
            F_within_active = F(active_start_idx:active_end_idx);
            F_within_active_abs = abs(F_within_active);
            maxF = max(F_within_active_abs);
        end

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

            P = c * v0 * exp(-(Ea / (kB * T))); % transition probability

            vD = d * P * (exp((z * e * d * F) / (2 * kB * T)) ...
                    - exp(-(z * e * d * F) / (2 * kB * T)));
            vD0 = d * P * ((z * e * d * F)/(kB * T));

            mu = vD ./ F; % mobility (drfit velocity divide field potential)
            mu0 = vD0 ./ F; % mobility at low electric field (dF/kBT << 1)
            mun = mu ./ mu0;

            % The non-linear transport 
            % non-linear --> longer retention time
            % just try different range, verigy w/ experimental staff
            % Plot: mun, Ea, thickness
        end
    end
end

