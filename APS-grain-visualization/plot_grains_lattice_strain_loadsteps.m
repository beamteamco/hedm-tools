% Plot the lattice strain in grain over all loading steps
% Loops over all grains ID'ed @ zero load
% Load all grain data
load APS_2016_s1_cycle_1_reconstruction_data;
plot_type = 4;

if(plot_type == 1)
    % Loop over all grains
    grain_loading_steps = [1:8]';
    grain_lat_strain = zeros(8, 1);
    hold on
    h1 = plot(grain_loading_steps, grain_lat_strain, 'LineWidth', 2, 'Color', 'b');
    set(gca, 'FontSize', 16);
    xlabel('Loading step');
    ylabel('Lattice strain along loading axis');
    ylim([-0.3 0.6]);
    xlim([1 8]);
    h1.XDataSource = 'grain_loading_steps';
    h1.YDataSource = 'grain_lat_strain';
    %
    for ii = 1:size(grain_correspondence_table, 1)
        grain_corr_id = grain_correspondence_table(ii, :);
        grain_lat_strain = zeros(numel(grain_corr_id), 1);
        for jj = 1:numel(grain_corr_id)
            if(grain_correspondence_table(ii, jj) ~= 0)
                grain_lat_strain(jj) = (APS_2016_s1_cycle_1_lattice_strain_grains{jj}(grain_correspondence_table(ii, jj), 2))/1.0e4;
            else
                grain_lat_strain(jj) = NaN;
            end
        end
        % Plot
        h2 = text(1.2, 0.55, ['Parent trans strain (predicted) = ' num2str(100*APS_2016_s1_cycle_1_grain_neighborhood_data(ii, 4), 4)]);
        h3 = text(1.2, 0.5,  ['Neighbors trans strain (predicted) = ' num2str(100*APS_2016_s1_cycle_1_grain_neighborhood_data(ii, 6), 4)]);
        h4 = text(1.2, 0.45, ['Number of neighbors = ' num2str(APS_2016_s1_cycle_1_grain_neighborhood_data(ii, 5))]);
        h5 = text(5, 0.55, ['Radial distance = ' num2str(APS_2016_s1_cycle_1_grain_neighborhood_data(ii, 3), 4)]);
        h6 = text(5, 0.5, ['Radius = ' num2str(APS_2016_s1_cycle_1_grain_neighborhood_data(ii, 2), 4)]);
        refreshdata
        print(gcf, ['grain_lat_strain_plots/grain_lat_strain_' num2str(ii)], '-dpng');
        delete([h2 h3 h4 h5 h6])
    end
    hold off;
elseif(plot_type == 2)
    % Plot lattice strains on a single plot
    to_delete = zeros(652, 1);
    for ii = 1:652
        if(min(grain_correspondence_table(ii, :)) == 0)
            to_delete(ii) = 1;
        end
    end
    grain_correspondence_table(to_delete == 1, :) = [];
    %
    grain_loading_steps = [1:8]';
    grain_lat_strain = NaN*zeros(8, 1);
    colors = distinguishable_colors(652);
    h = plot(grain_loading_steps, grain_lat_strain, 'LineWidth', 2, 'Color', 'b');
    set(gca, 'FontSize', 16);
    xlabel('Loading step');
    ylabel('Lattice strain along loading axis');
    ylim([-0.3 0.6]);
    xlim([1 8]);
    hold on;
    set(gcf, 'Visible', 'off');
    %
    for ii = 1:size(grain_correspondence_table, 1)
        grain_corr_id = grain_correspondence_table(ii, :);
        grain_lat_strain = zeros(numel(grain_corr_id), 1);
        for jj = 1:numel(grain_corr_id)
            if(grain_correspondence_table(ii, jj) ~= 0)
                grain_lat_strain(jj) = (APS_2016_s1_cycle_1_lattice_strain_grains{jj}(grain_correspondence_table(ii, jj), 2))/1.0e4;
            else
                grain_lat_strain(jj) = NaN;
            end
        end
        plot(grain_loading_steps, grain_lat_strain, 'LineWidth', 1, 'Color', colors(ii, :));
    end
    set(gcf, 'Visible', 'on');
elseif(plot_type == 3)
    % Plot CM Y on a single plot
    to_delete = zeros(652, 1);
%     for ii = 1:652
%         if(min(grain_correspondence_table(ii, :)) == 0)
%             to_delete(ii) = 1;
%         end
%     end
    grain_correspondence_table(to_delete == 1, :) = [];
    %
    grain_loading_steps = [1:8]';
    grain_com_y = NaN*zeros(8, 1);
    colors = distinguishable_colors(652);
    h = plot(grain_loading_steps, grain_com_y, 'LineWidth', 2, 'Color', 'b');
    set(gca, 'FontSize', 16);
    xlabel('Loading step');
    ylabel('COM (Y)');
    %ylim([-0.3 0.6]);
    xlim([1 8]);
    hold on;
    set(gcf, 'Visible', 'off');
    %
    for ii = 1:size(grain_correspondence_table, 1)
        grain_corr_id = grain_correspondence_table(ii, :);
        grain_com_y = zeros(numel(grain_corr_id), 1);
        for jj = 1:numel(grain_corr_id)
            if(grain_correspondence_table(ii, jj) ~= 0)
                grain_com_y(jj) = (APS_2016_s1_cycle_1_com{jj}(grain_correspondence_table(ii, jj), 2));
            else
                grain_com_y(jj) = NaN;
            end
        end
        plot(grain_loading_steps, grain_com_y, 'LineWidth', 1, 'Color', colors(ii, :));
    end
    set(gcf, 'Visible', 'on');
elseif(plot_type == 4)
    % Plot strain triaxiality (strains(2)/sum(strains)) on a single plot
    to_delete = zeros(652, 1);
    for ii = 1:652
        if(min(grain_correspondence_table(ii, :)) == 0)
            to_delete(ii) = 1;
        end
    end
    grain_correspondence_table(to_delete == 1, :) = [];
    %
    grain_loading_steps = [1:8]';
    grain_triaxiality = NaN*zeros(8, 1);
    colors = distinguishable_colors(652);
    h = plot(grain_loading_steps, grain_triaxiality, 'LineWidth', 2, 'Color', 'b');
    set(gca, 'FontSize', 16);
    xlabel('Loading step');
    ylabel('Triaxiality');
    %ylim([-0.3 0.6]);
    xlim([1 8]);
    hold on;
    set(gcf, 'Visible', 'off');
    %
    for ii = 1:size(grain_correspondence_table, 1)
        grain_corr_id = grain_correspondence_table(ii, :);
        grain_triaxiality = zeros(numel(grain_corr_id), 1);
        for jj = 1:numel(grain_corr_id)
            if(grain_correspondence_table(ii, jj) ~= 0)
                grain_triaxiality(jj) = abs(APS_2016_s1_cycle_1_lattice_strain_grains{jj}(grain_correspondence_table(ii, jj), 2))/(abs(APS_2016_s1_cycle_1_lattice_strain_grains{jj}(grain_correspondence_table(ii, jj), 1)) + abs(APS_2016_s1_cycle_1_lattice_strain_grains{jj}(grain_correspondence_table(ii, jj), 2)) + abs(APS_2016_s1_cycle_1_lattice_strain_grains{jj}(grain_correspondence_table(ii, jj), 3)));
            else
                grain_triaxiality(jj) = NaN;
            end
        end
        plot(grain_loading_steps, grain_triaxiality, 'LineWidth', 1, 'Color', colors(ii, :));
    end
    set(gcf, 'Visible', 'on');
end