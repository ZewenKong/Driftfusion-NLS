function AM15 = lightsource(source_type, lambda)
    % loads light sources for use with BeerLambert function
    % lambda - range of wavelengths in nanometers

    if source_type == 'AM15'
        % Load AM1.5
        AM15_data = readtable('AM15.xls', 'VariableNamingRule', 'preserve');
        AM15 = 1e-3 * interp1(AM15_data.(1), AM15_data.(2), lambda, 'linear', 'extrap'); % Converts to W cm-2 nm-1
    end

end
