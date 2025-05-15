classdef dfcalc

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

            j_bv = j0 * (exp((1 - alpha_e) * F * (E - E_eq)/ (R * T)) - exp((- alpha_e) * F * (E - E_eq)/ (R * T)));
        end

        function save_calcdata(E_boundary, E_eq, eta, a_boundary, t)
            persistent calcdata;

            if isempty(calcdata)
                calcdata.E_boundary = [];
                calcdata.E_eq = [];
                calcdata.eta = [];
                calcdata.a_boundary = [];
                calcdata.t = [];
            end
            
            calcdata.E_boundary(end + 1, 1) = E_boundary;
            calcdata.E_eq(end + 1, 1) = E_eq;
            calcdata.eta(end + 1, 1) = eta;
            calcdata.a_boundary(end + 1, 1) = a_boundary;
            calcdata.t(end + 1, 1) = t;

            assignin('base', 'calcdata', calcdata);
        end
    end
end