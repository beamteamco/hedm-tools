% Visualize individual Grains.csv file
grains_file_name = '/Volumes/Secondary HD/Work/Research_data/HEDM_Study_of_Deformation_in_SMA/CHESS_Dec_2015/results/MIDAS/Grains.csv';                 % grains file name
make_plots = true;
% Read Grains.csv data
grains = load(grains_file_name);

% Store useful quantities in arrays
% Note: after applying the RESRF2APS transformation, loading axis is Y,
% just like heXRD
RESRF2APS   = RMatOfQuat(QuatOfESRF2APS);
%
num_grains     = size(grains, 1);
%
grains_com = RESRF2APS*grains(:, 11:13)';
grains_com = grains_com';
grains_lattice_params = grains(:, 14:19);
grains_strains_fab = grains(:, 25:33);
% After applying ESRF (FABLE) -> APS transformation, grain strain
% tensor in the laboratory coordinates. Y (2, 2) is loading axis.
for i = 1:1:num_grains
    grains_strains_fab_matrix(:,:,i)   =  RESRF2APS*reshape(grains_strains_fab(i, :), 3, 3)'*RESRF2APS';
end
grains_strains_ken = grains(:, 34:42);
grains_orient = grains(:, 2:10);
grains_radius = grains(:, 23);
grains_completeness = grains(:, 24);
grains_radial_dist = sqrt(grains_com(:, 1).^2 + grains_com(:, 3).^2);

disp(' ')
disp(['Total number of grains = ' num2str(size(grains, 1))])
disp(['Mean radius = ' num2str(mean(grains_radius, 1))])
disp(['Mean strains (FABLE XYZ, Y = loading) = ' num2str(mean(squeeze(grains_strains_fab_matrix(1, 1, :)))) ', ' num2str(mean(squeeze(grains_strains_fab_matrix(2, 2, :)))) ', ' num2str(mean(squeeze(grains_strains_fab_matrix(3, 3, :))))])
disp(['standard dev strains (FABLE XYZ, Y = loading) = ' num2str(std(squeeze(grains_strains_fab_matrix(1, 1, :)))) ', ' num2str(std(squeeze(grains_strains_fab_matrix(2, 2, :)))) ', ' num2str(std(squeeze(grains_strains_fab_matrix(3, 3, :))))])
disp(' ')

% Calculate the orientation (MTEX) for each grain
for i = 1:1:num_grains
    RMats(:,:,i)   =  RESRF2APS*reshape(grains(i,2:10), 3, 3)';
end
qsym    = CubSymmetries; Rsym    = RMatOfQuat(qsym);
quat    = ToFundamentalRegionQ(QuatOfRMat(RMats), qsym);

cs = crystalSymmetry('cubic');
ss = specimenSymmetry('triclinic');
q = quaternion(quat);
r = rotation(q);
o = orientation(q, cs, ss);
%
oM = ipdfHSVOrientationMapping(crystalSymmetry('432'));
% Colors are generated for [0 0 1] IPF. Need to fix that.
oM.inversePoleFigureDirection = vector3d(0,1,0);
rgb = oM.orientation2color(o);
% Rotate stiffness tensor to global coords in each gran
% Calculate stiffness along loading direction (in MPa)
%     for i = 1:num_grains
%         C_rot = rotate(stiffness_tensor, r(i));
%         stiffness_loading(i) = YoungsModulus(C_rot , yvector)/1e6;
%     end
if(make_plots)
    % Visualize IPF
    plotIPDF(o, yvector, 'Property', squeeze(grains_strains_fab_matrix(2, 2, :)), 'xAxisDirection','east', 'MarkerSize', 4);
    colorbar; colormap jet;
    % Visualize COM
    figure;
    scatter3(grains_com(:, 1), grains_com(:, 2), grains_com(:, 3), 121, grains_radial_dist > 400, 'filled');
    campos([0 1800 0]); camtarget([0 200 0]); camup([0 0 1]);
    axis equal;
    axis vis3d;
end
%