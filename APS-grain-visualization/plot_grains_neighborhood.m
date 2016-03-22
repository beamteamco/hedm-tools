% Plot 3D tesselation at each loading step by tracking "vanishing" grains
% Load all grain data
load APS_2016_s1_cycle_1_reconstruction_data;
% Neighborhood data
% id radius radial_dist tens_schmid_plate num_of_neighbors avg_tens_schmid_plate_neighbors
APS_2016_s1_cycle_1_grain_neighborhood_data = zeros(652, 6);
for ii = 1:652
    parent_grain_id = ii;  % Neighbors of this grain will be plotted
    neighbor_distance_threshold = 150; % Set to 1 to plot just parent, ~ 100 for nearest neighbors
    make_plot = 0;
    % Make a copy of key grain/tess parameters
    elements_ii = elements;
    grain_color_ii = grain_color;
    grain_id_ii = grain_id;
    grains_com_ii = APS_2016_s1_cycle_1_com{1};
    elements_to_remove = zeros(size(elements, 1), 1);
    % Loop over grains
    parent_com = grains_com_ii(parent_grain_id, :);
    for jj = 1:size(grains_com_ii, 1)
        % FInd dist between parent and current grain from COM data
        com_ii = grains_com_ii(jj, :);
        dist_parent_jj = sqrt(sum((com_ii - parent_com).*(com_ii - parent_com)));
        % If dist > threshold then delete that grain
        if(dist_parent_jj > neighbor_distance_threshold)
            elements_to_remove(grain_id(:, 1) == jj) = 1;
        end
    end
    % Actually remove elements etc for distant grains
    elements_ii(elements_to_remove == 1, :) = [];
    grain_color_ii(elements_to_remove == 1, :) = [];
    grain_id_ii(elements_to_remove == 1, :) = [];
    
    % Calculate neighborhood parameters
    APS_2016_s1_cycle_1_grain_neighborhood_data(ii, 1) = ii;
    APS_2016_s1_cycle_1_grain_neighborhood_data(ii, 2) = APS_2016_s1_cycle_1_radius_grains{1}(ii);
    APS_2016_s1_cycle_1_grain_neighborhood_data(ii, 3) = sqrt(grains_com_ii(ii, 1)*grains_com_ii(ii, 1) + grains_com_ii(ii, 3)*grains_com_ii(ii, 3));
    APS_2016_s1_cycle_1_grain_neighborhood_data(ii, 4) = TPSCHMID(ii);
    APS_2016_s1_cycle_1_grain_neighborhood_data(ii, 5) = numel(unique(grain_id_ii));
    APS_2016_s1_cycle_1_grain_neighborhood_data(ii, 6) = mean(TPSCHMID(unique(grain_id_ii)));
    % Plot
    if(make_plot)
        figure;
        scatter3(elements_ii(:, 1), elements_ii(:, 2), elements_ii(:, 3), 64, grain_color_ii, 'filled', 'marker', 's')
        set(gca, 'FontSize', 16);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        xlim([-0.55 0.55]); ylim([-0.55 0.55]); zlim([-0.05 0.75]);
        axis vis3d; grid off; box on;
    end
end

