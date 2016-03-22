% Plot 3D tesselation at each loading step by tracking "vanishing" grains
% Load all grain data
load APS_2016_s1_cycle_1_reconstruction_data;
% Loop over all loading steps
for ii = 1:numel(APS_2016_s1_cycle_1_com)
    elements_ii = elements;
    grain_color_ii = grain_color;
    elements_to_remove = zeros(size(elements, 1), 1);
    for jj = 1:size(grain_correspondence_table, 1)
        if(grain_correspondence_table(jj, ii) == 0)
            elements_to_remove(grain_id(:, 1) == jj) = 1;
        end
    end
    elements_ii(elements_to_remove == 1, :) = [];
    grain_color_ii(elements_to_remove == 1, :) = [];
    %
    figure;
    scatter3(elements_ii(:, 1), elements_ii(:, 2), elements_ii(:, 3), 64, grain_color_ii, 'filled', 'marker', 's')
    title(['Loading step ' num2str(ii)]);
    set(gca, 'FontSize', 16);
    xlabel('X'); ylabel('Y'); zlabel('Z'); 
    xlim([-0.55 0.55]); ylim([-0.55 0.55]); zlim([-0.05 0.75]);
    print(gcf, ['tess_loading_step_' num2str(ii)], '-dpng');
    close all;
end