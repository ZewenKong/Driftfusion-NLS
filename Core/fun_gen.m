function fun = fun_gen(fun_type)
% 
% Function generator for generation profile
% Some of the functions in this code are based on SFG.M by Hiroyuki Kato Copyright (c) 2013
% The license details are contained with sfg_license
% G_FUN_TYPE = Type of function
%
% 'constant'
% COEFF = [Amplitude]
%
% 'sweep'
% COEFF = [Amplitude_initial, Amplitude_final, tmax]
% Here A is ignored - this was the easiest way to maintain backwards
% compatibility with different protocols and may be updated in future
% releases
%
% 'square'
% COEFF = [A_low, A_high, time_period, duty_cycle]
%
% 'sin'
% COEFF = [DC_offset, Delta_AC, frequency, phase]
%
% 'sweepAndStill'
% COEFF = [A_start, A_end, sweep_duration]
%
% LICENSE
% Copyright (C) 2020  Philip Calado, Ilario Gelmetti, and Piers R. F. Barnes
% Imperial College London
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
%% - - - - - - - - - - CODE START - - - - - - - - - -

switch fun_type
    case 'constant'
        fun = @(coeff, t) coeff(1); 
        % constant voltage value through simulation
        % coeff(1), voltage value
        %
    case 'sweep'
        fun = @(coeff, t) coeff(1) + (coeff(2) - coeff(1)) * t / coeff(3);
        % coeff(1), start voltage
        % coeff(2), end voltage
        % coeff(3), sweeping time
        %
    case 'square'
        fun = @(coeff, t) coeff(1) + (coeff(2) - coeff(1)) * lt(mod(t, coeff(3)) * 1/coeff(3), coeff(4)/100);
        % square pulse (generate intensity array)
        % coeff(1), low voltage (baseline)
        % coeff(2), high voltage
        % coeff(3), square pulse period
        % coeff(4), duty cycle in percent
        %
    case 'sin'
        fun = @(coeff, t) coeff(1) + coeff(2) * (sin(2 * pi * coeff(3) * t + coeff(4)));
        % coeff(1), DC offset
        % coeff(2), amplitude
        % coeff(3), frequency in Hz
        % coeff(4), phase shift in radians
        %
    case 'tri'
        fun = @(coeff, t) triangle_fun(coeff, t);
        % COEFF = [OFFSET, V1, V2, periods, tperiod]  
        % tmax is defined by the input
        %
    case 'smoothed_square'
        fun = @(coeff, t) smoothed_square(coeff, t);
        % COEFF = [A_start, A_pulse, period, duty_cycle]  
        %
    case 'sweepAndStill'
        fun = @(coeff, t) coeff(2) + (coeff(1)-coeff(2))*max(0,1-t/coeff(3));
        % COEFF = [A_start, A_end, sweep_duration]
        % do a sweep from coeff(1) to coeff(2) with the duration of
        % coeff(3), then stay stable at coeff(2)
end