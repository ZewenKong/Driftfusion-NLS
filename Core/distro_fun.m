classdef distro_fun

    % DRIFTFUSION distribution function class calculates carrier densities for
    % different distribution functions.

    %% - - - - - - - - - - CODE START - - - - - - - - - -

    properties (Constant)

        kB = 8.617330350e-5; % Boltzmann constant [eV K^-1]
        uplimit = 10; % Upper limit for the integration [eV]
    end

    methods (Static)

        function n = nfun(Nc, Ec, Efn, par)

            % nfun, negative function (electron density, conduction band
            % energy, electron quasi fermi level, par)
            % n, charge carrier density

            kT = par.kB * par.T;

            switch par.prob_distro_function

                case 'Fermi'

                    % Fermi dirac integral for obtaining electron densities
                    % Nc = conduction band density of states
                    % Ec = conduction band energy
                    % Ef = Fermi level
                    % T = temperature
                    % See Schubert 2015, pp. 130

                    n = zeros(1, length(Nc));

                    for i = 1:length(Nc)

                        if isnan(Nc(i)) == 0 % ignores interfaces

                            fn = @(E) ((E / kT) .^ 0.5) ./ (1 + exp((E - Efn(i) + Ec(i)) / kT));
                            n(i) = real(((2 * Nc(i)) / (kT * pi ^ 0.5)) * integral(fn, 0, distro_fun.uplimit));
                        end

                    end

                case 'Blakemore'

                    eta_n = (Efn - Ec) ./ kT;
                    n = Nc .* (1 ./ (exp(-eta_n) + par.gamma));

                case 'Boltz' % Boltzmann approximation is used in the this study (RPP based memristor)

                    % Nc, conduction band density of states (comes from effective mass and Planck's const)
                    % Efn, electron quasi fermi level
                    % Ec, conduction band energy

                    n = Nc .* exp((Efn - Ec) ./ kT);
            end

        end

        function p = pfun(Nv, Ev, Efp, par)

            kT = par.kB * par.T;

            switch par.prob_distro_function

                case 'Fermi'

                    % Fermi dirac integral for obtaining electron densities
                    % Nc = conduction band density of states
                    % Ec = conduction band energy
                    % Ef = Fermi level
                    % T = temperature
                    % see Schubert 2015, pp. 130

                    p = zeros(1, length(Nv));

                    % reflecting the energy makes the integral easier for some
                    % reason- doesn't seem to like integrating from negative infinitiy...

                    Efp = Ev - (Efp - Ev);

                    for i = 1:length(Nv)

                        if isnan(Nv(i)) == 0 % ignores interfaces
                            fp = @(E) ((E / kT) .^ 0.5) ./ (1 + exp((E - Efp(i) + Ev(i)) / kT));
                            p(i) = real(((2 * Nv(i)) / (kT * pi ^ 0.5)) * integral(fp, 0, distro_fun.uplimit));
                        end

                    end

                case 'Blakemore'

                    eta_p = (Ev - Efp) ./ kT;
                    p = Nv .* (1 ./ (exp(-eta_p) + par.gamma));

                case 'Boltz'

                    p = Nv .* exp((Ev - Efp) ./ kT);
            end

        end

        function Dnfd = Dn_fd_fun(Nc, Ec, Efn, mu_n, T)

            if isnan(Ec) == 0 % ignores interfaces

                % calculates Fermi Dirac diffusion coefficient function
                %
                % currently uses upper limit for the integral of CB +16eV- using
                % higher values causes problems resulting in Nan results for Dn
                % - requires further investigation
                % Nc = conduction band density of states
                % Ec = conduction band energy
                % Ef = Fermi level
                % T = temperature
                % see Schubert 2015, pp. 130
                kB = 8.617330350e-5; % Boltzmann constant [eV K^-1]
                kT = kB * T;
                e = 1.61917e-19; % elementary charge in Coulombs.

                for i = 1:length(Efn)

                    f = @(E) 1 ./ (1 + exp((E - Efn(i) + Ec) / kT)); % Fermi- Dirac function
                    dfdE = @(E) exp((E - Efn(i) + Ec) / kT) ./ (kT * (exp((E - Efn(i) + Ec) / kT) + 1) .^ 2);

                    g = @(E) (E / kT) .^ 0.5; % DOS function based on 3D semiconductor
                    h = @(E) g(E) .* f(E);
                    k = @(E) g(E) .* dfdE(E);

                    n(i) = real(((2 * Nc) / (kT * pi ^ 0.5)) * integral(f, 0, distro_fun.uplimit));
                    dndE(i) = ((2 * Nc) / (kT * pi ^ 0.5)) * integral(dfdE, 0, distro_fun.uplimit);

                end

                Dnfd.Dnfun = mu_n * (n ./ dndE);
                Dnfd.n_fd = n;
                Dnfd.Efn = Efn;

            end

        end

        function Dpfd = Dp_fd_fun(Nv, Ev, Efp, mu_p, T)

            if isnan(Ev) == 0 % ignores interfaces

                % calculates Fermi Dirac diffusion coefficient function
                % currently uses upper limit for the integral of CB +16eV- using
                % higher values causes problems resulting in Nan results for Dn
                % - requires further investigation
                % Nc = conduction band density of states
                % Ec = conduction band energy
                % Ef = Fermi level
                % T = temperature
                % see Schubert 2015, pp. 130

                % reflecting the energy makes the integral easier for some
                % reason- doesn't seem to liek integrating from negative infinitiy...
                Efp_trans = 2 * Ev - Efp;

                kB = 8.617330350e-5; % Boltzmann constant [eV K^-1]
                kT = kB * T;
                e = 1.61917e-19; % elementary charge in Coulombs.

                for i = 1:length(Efp)
                    f = @(E) 1 ./ (1 + exp((E - Efp_trans(i) + Ev) / kT)); % Fermi- Dirac function
                    dfdE = @(E) exp((E - Efp_trans(i) + Ev) / kT) ./ (kT * (exp((E - Efp_trans(i) + Ev) / kT) + 1) .^ 2);

                    g = @(E) (E / kT) .^ 0.5; % DOS function based on 3D semiconductor
                    h = @(E) g(E) .* f(E);
                    k = @(E) g(E) .* dfdE(E);

                    p(i) = real(((2 * Nv) / (kT * pi ^ 0.5)) * integral(f, 0, distro_fun.uplimit));
                    dpdE(i) = ((2 * Nv) / (kT * pi ^ 0.5)) * integral(dfdE, 0, distro_fun.uplimit);
                end

                Dpfd.Dpfun = mu_p * (p ./ dpdE);
                Dpfd.p_fd = p;
                Dpfd.Efp = Efp;

            end

        end

        function Dsol = Dnlook(n, Dnfun, n_fd)

            if n > max(n_fd)
                warning('Carrier density is outside of specified limits for the Fermi Dirac integral: Try increasing the FERMI_LIMIT in PC')
            end

            % to avoid problems with negative carrier densities
            if n < 0
                n = 0;
            end

            if n > min(n_fd)
                pp = find(n_fd <= n);
                pp = pp(end);
                Dsol = Dnfun(1, pp);
                % Dsol = interp1(n_fd, Dnfun(1,pp), n, 'pchip');%'nearest', 'extrap');
            else
                Dsol = min(Dnfun);
            end

            %
        end

        function Dsol = Dplook(p, Dpfun, p_fd)

            if p > max(p_fd)
                warning('Carrier density is outside of specified limits for the Fermi Dirac integral: Try increasing the FERMI_LIMIT in PC')
            end

            % to avoid problems with negative carrier densities
            if p < 0
                p = 0;
            end

            if p > min(p_fd)
                pp = find(p_fd >= p);
                pp = pp(end);
                Dsol = Dpfun(1, pp);
                % Dsol = interp1(p_fd, Dpfun(1,pp), p, 'pchip');
            else
                Dsol = min(Dpfun);
            end

        end

        function Efn_fd = Efn_fd_fun(n, Efn, n_fd)

            if n > min(n_fd)
                pp = find(n_fd <= n);
                pp = pp(end);
                Efn_fd = Efn(1, pp);
            else
                Efn_fd = min(Efn);
            end

        end

        function Efp_fd = Efp_fd_fun(p, Efp, p_fd)

            if p > min(p_fd)
                pp = find(p_fd <= p);
                pp = pp(1);
                Efp_fd = Efp(1, pp);
            else
                Efp_fd = min(Efp);
            end

        end

    end

end
