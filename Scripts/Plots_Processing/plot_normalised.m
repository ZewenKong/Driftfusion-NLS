h = openfig('Model_Validation/SSP_Model/ana_results_mun/lat.fig');

% Set figure size to be square
size_inch = 5; % for example, 6 inches x 6 inches

set(h, 'Units', 'inches', 'Position', [1 1 size_inch size_inch]);

set(h, 'PaperPositionMode', 'auto');
print(h, 'your_file.png', '-dpng', '-r300');