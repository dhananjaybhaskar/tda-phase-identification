% Author: Dhananjay Bhaskar
% Plot heatmap of particle density at varying adhesion and self-propulsion force

num_neighbor_mat = zeros(11, 11);

adhesion_values = linspace(0.05, 0.25, 11);
polarity_values = linspace(0.005, 0.025, 11);

xlabelvec = {};
ylabelvec = {};

for i = 1 : length(adhesion_values)
    xlabelvec{end+1} = num2str(adhesion_values(i));
end

for i = 1 : length(polarity_values)
    ylabelvec{end+1} = num2str(polarity_values(i));
end

addpath(strcat(pwd, filesep, 'cbrewer'))
CT = cbrewer('seq', 'Reds', 121, 'PCHIP');

for adh_idx = 1 : 11

    for pol_idx = 1 : 11

        folder_name = strcat(num2str(adh_idx), '_adh_', num2str(pol_idx), '_pol');

        nbdfile = 'Neighbors_200000.dat';

        if ~exist(strcat(folder_name, filesep, nbdfile), 'file')
            continue;
        end

        disp(folder_name)

        nbd = csvread(strcat(folder_name, filesep, nbdfile));
        numcells = size(nbd, 1);

        frame_avg_num_neighbors = 0;
        neighbors = [];

        for j = 1 : numcells

            frame_avg_num_neighbors = frame_avg_num_neighbors + nbd(j);
            neighbors = [neighbors nbd(j)];

        end

        frame_avg_num_neighbors = frame_avg_num_neighbors/numcells;
        average_neighbors = mean(neighbors);

        disp(strcat('Adhesion:', num2str(adh_idx), ' Polarity:', num2str(pol_idx), ' # Neighbors:', num2str(average_neighbors)));

        num_neighbor_mat(adh_idx, pol_idx) = average_neighbors;

    end

end

save('local_density_noprolif.mat', 'num_neighbor_mat')

h = heatmap(flipud(num_neighbor_mat'), 'colormap', CT, 'GridVisible', 'off', 'CellLabelColor', 'none');

h.XLabel = 'Adhesion Parameter';
h.YLabel = 'Random Self-Propulsion Force';
h.YDisplayLabels = flip(ylabelvec);
h.XDisplayLabels = xlabelvec;

set(h, 'FontSize', 11, 'FontName', 'Myriad Pro Regular')

set(gca, 'XColor', 'black');
set(gca, 'YColor', 'black')
set(h, 'LineWidth', 1);

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 6.25 6.25]);

print('num_neighbors_heatmap.eps', '-depsc')
