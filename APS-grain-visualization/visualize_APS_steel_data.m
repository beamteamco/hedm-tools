% APS Steel data visualization
% Columns
% Sp_ID O[0][0] O[0][1] O[0][2] O[1][0] O[1][1] O[1][2] O[2][0] O[2][1] O[2][2] X Y Z a b c alpha beta gamma Err1 Err2 Err3 MeanRadius Confidence 
%
grains  = importdata('~/Desktop/Grains.csv');
search_volume = 2 * 2 * 2;
sample_volume = 1 * 0.1 * 0.5;

grains = grains.data;

RESRF2APS   = RMatOfQuat(QuatOfESRF2APS);

nGrains     = size(grains, 1);
for i = 1:1:nGrains
    RMats(:,:,i)   =  RESRF2APS*reshape(grains(i,2:10), 3, 3)';
end
qsym    = CubSymmetries; Rsym    = RMatOfQuat(qsym);
quat    = ToFundamentalRegionQ(QuatOfRMat(RMats), qsym);
rod     = RodOfQuat(quat);

grains_COM = grains(:, [11 12 13]);
grains_COM = (RESRF2APS * grains_COM')';

grains_radius = grains(:, 23);
grains_confidence = grains(:, 24);
grains_radius_scale_factor = sample_volume / search_volume;
grains_radius_scaled = grains_radius * grains_radius_scale_factor;

scatter3(grains_COM(:, 1), grains_COM(:, 2), grains_COM(:, 3), grains_radius, grains_confidence, 'filled');
%xlim([-500 500]);
%ylim([-250 250]);
%zlim([-50 50]);
xlabel('X'); ylabel('Y'); zlabel('Z');

axis vis3d; axis equal; grid off;