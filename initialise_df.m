function initialise_df

    % plotting defaults
    set(0, 'DefaultLineLinewidth', 1.5);
    set(0, 'DefaultAxesFontSize', 16);
    set(0, 'DefaultFigurePosition', [600, 400, 600, 450]);
    set(0, 'DefaultAxesXcolor', [0, 0, 0]);
    set(0, 'DefaultAxesYcolor', [0, 0, 0]);
    set(0, 'DefaultAxesZcolor', [0, 0, 0]);
    set(0, 'DefaultTextColor', [0, 0, 0]);

    % add file paths to functions, avoid hidden directories like .git
    addpath(pwd)
    files = dir(pwd);
    dirs = files(cat(1, files.isdir));

    for i = 1:length(dirs)
        dirname = dirs(i).name;

        if ~any(regexp(dirname, '^\.'))
            addpath(dirname);
        end

    end

end
