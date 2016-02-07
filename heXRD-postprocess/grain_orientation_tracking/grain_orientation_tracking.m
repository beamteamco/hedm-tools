% Reads a list of 'accepted_orientations.dat' files and determines the
% correspondence between a single grain at different loading steps
%
% List of heXRD orientation files to read
accepted_orientations_list = {'accepted_orientations_18.dat', ...
                              'accepted_orientations_25.dat'};
track_grains_com = true;
grains_list = {'grains_18.out', ...
              'grains_25.out'};
% Crystal symmetry (cubic, monoclinic etc) of the material (This sent to MTEX)
crystal_symmetry = 'cubic';
% Grains with misorientaion less than the threshold are considered same
misorientation_threshold = 1;
% Read orientations for the reference grain list
a_root = importdata(accepted_orientations_list{1});
q_root = quaternion(a_root');
o_root = orientation(q_root, crystalSymmetry(crystal_symmetry), specimenSymmetry('triclinic'));
%
if(track_grains_com)
    grains_root = importdata(grains_list{1});
    grains_root = grains_root.data;
end
% Create tables to save grain correspondence, misorientation for the
% corresponding grains, COM difference etc
grain_correspondence_table = zeros(size(a_root, 1), numel(accepted_orientations_list));
grain_correspondence_table(:, 1) = [1:size(a_root, 1)]';
grain_misorientation_table = zeros(size(a_root, 1), numel(accepted_orientations_list));
grain_misorientation_table(:, 1) = [1:size(a_root, 1)]';
grain_distance_table = zeros(size(a_root, 1), numel(accepted_orientations_list));
grain_distance_table(:, 1) = [1:size(a_root, 1)]';
% Loop over the rest of the orientation files
for ii = 2:numel(accepted_orientations_list)
    % Read orientation info for the current list
    a_current = importdata(accepted_orientations_list{ii});
    q_current = quaternion(a_current');
    o_current = orientation(q_current, crystalSymmetry(crystal_symmetry), specimenSymmetry('triclinic'));
    %
    if(track_grains_com)
        grains_current = importdata(grains_list{ii});
        grains_current = grains_current.data;
    end
    % Calculate the pairwise misorientation between grains from the two
    % lists
    grain_misorientation = pdist2(reshape(o_root, size(o_root, 2), 1), reshape(o_current, size(o_current, 2), 1), @angle) * 180.0 / pi;
    if(track_grains_com)
        grain_distance = pdist2(grains_root(:, [7 8 9]), grains_current(:, [7 8 9]));
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
    % Save the corresponding grain list
    grain_correspondence_table(:, ii) = min_misorientation_index;
    grain_misorientation_table(:, ii) = min_misorientation;
    if(track_grains_com)
        grain_distance_table(:, ii) = min_distance;
    end
end