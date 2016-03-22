% Reads a list of 'Grains.csv' files and determines the
% correspondence between a single grain at different loading steps
%
% List of MIDAS grains files to read
track_grains_com = true;
orientation_mat_file = 'APS_2016_s1_cycle_1_reconstruction_data';
y_correction = [0 0 15 30 45 30 15 -15];
% Crystal symmetry (cubic, monoclinic etc) of the material (This sent to MTEX)
crystal_symmetry = 'cubic';
% Grains with misorientaion less than the threshold are considered same
misorientation_threshold = 1;
distance_threshold = 70;
% Read orientations
load(orientation_mat_file);
o_root = APS_2016_s1_cycle_1_orientations{1};
c_root = APS_2016_s1_cycle_1_com{1}; % center of mass
%
% Create tables to save grain correspondence, misorientation for the
% corresponding grains, COM difference etc
grain_correspondence_table = zeros(size(o_root, 2), numel(APS_2016_s1_cycle_1_orientations));
grain_correspondence_table(:, 1) = [1:size(o_root, 2)]';
grain_misorientation_table = zeros(size(o_root, 2), numel(APS_2016_s1_cycle_1_orientations));
grain_misorientation_table(:, 1) = [1:size(o_root, 2)]';
grain_distance_table = zeros(size(o_root, 2), numel(APS_2016_s1_cycle_1_orientations));
grain_distance_table(:, 1) = [1:size(o_root, 2)]';
% Loop over the rest of the orientation files
for ii = 2:numel(APS_2016_s1_cycle_1_orientations)
    % Read orientation info for the current list
    o_current = APS_2016_s1_cycle_1_orientations{ii};
    c_current = APS_2016_s1_cycle_1_com{ii};  % center of mass
    c_current = [c_current(:, 1) c_current(:, 2) - y_correction(ii) c_current(:, 3)];
    %
    % Calculate the pairwise misorientation between grains from the two
    % lists
    grain_misorientation = pdist2(reshape(o_root, size(o_root, 2), 1), reshape(o_current, size(o_current, 2), 1), @angle) * 180.0 / pi;
    if(track_grains_com)
        grain_distance = pdist2(c_root, c_current);
    end
    [min_misorientation, min_misorientation_index] = min(grain_misorientation, [], 2);
    if(track_grains_com)
        for jj = 1:numel(min_misorientation_index)
            min_distance(jj) = grain_distance(jj, min_misorientation_index(jj));
        end
    end
    % Remove misorientations above the threshold
    min_misorientation_index(min_misorientation > misorientation_threshold) = 0;
    if(track_grains_com)
        min_distance(min_misorientation > misorientation_threshold) = NaN;
    end
    min_misorientation(min_misorientation > misorientation_threshold) = NaN;
    % Remove misorientations above distance threshold
    if(track_grains_com)
        min_misorientation_index(min_distance > distance_threshold) = 0;
        min_misorientation(min_distance > distance_threshold) = NaN;
        if(track_grains_com)
            min_distance(min_distance > distance_threshold) = NaN;
        end
    end
    % Save the corresponding grain list
    grain_correspondence_table(:, ii) = min_misorientation_index;
    grain_misorientation_table(:, ii) = min_misorientation;
    if(track_grains_com)
        grain_distance_table(:, ii) = min_distance;
    end
end