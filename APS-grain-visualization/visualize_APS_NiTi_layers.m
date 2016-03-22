% Visualize MIDAS grain reconstruction
%    Read Grains.csv file for multiple layers and plot COM, orientations
%

APS_2016_s1_cycle_1_radius_grains = {};
APS_2016_s1_cycle_1_lattice_strain = [];
APS_2016_s1_cycle_1_lattice_strain_surf = [];
APS_2016_s1_cycle_1_lattice_strain_stdev = [];
APS_2016_s1_cycle_1_lattice_strain_stdev_surf = [];
APS_2016_s1_cycle_1_num_of_grains = [];
APS_2016_s1_cycle_1_num_of_grains_surf = [];
APS_2016_s1_cycle_1_orientations = {};
APS_2016_s1_cycle_1_com = {};
APS_2016_s1_cycle_1_lattice_strain_grains = {};

for ll = [0:6 8]
    % INPUTS
    % Load stiffness tensor for NiTi from file
    stiffness_tensor = loadTensor('NiTi_cubic_elastic_constants.data' , crystalSymmetry('cubic'), 'name', 'NiTi stiffness');
    total_num_layers = 5;
    layer_thickness  = 150;
    grains_file_stem = ['grains_files/Grains_NiTi_s1_0' num2str(ll) '_Layer'];                 % grains file name = stem + layer_num
    grains_to_plot = [];                                                   % Array to plot specific grains, empty to plot all
    make_plots = false;
    % Read Grains.csv data
    grains = [];
    grains_layer_id = [];
    for layer_num = 1:total_num_layers
        grains_layer = load([grains_file_stem num2str(layer_num) '.csv']);
        grains_layer(:, 13) = grains_layer(:, 13) + layer_thickness*(layer_num - 1);
        grains = [grains; grains_layer];
        grains_layer_id = [grains_layer_id; layer_num * ones(size(grains_layer, 1), 1)];
    end
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
    disp(['Loading step ' num2str(ll)])
    disp(['Total number of grains = ' num2str(size(grains, 1))])
    disp(['Grains with radial dist > 300 = ' num2str(sum(grains_radial_dist > 300))])
    disp(['Fraction of surface grains = ' num2str(sum(grains_radial_dist > 300)/size(grains, 1))])
    disp(['Mean radius = ' num2str(mean(grains_radius, 1))])
    disp(['Mean radius (surf) = ' num2str(mean(grains_radius(grains_radial_dist(:) > 300), 1))])
    disp(['Mean strains (FABLE XYZ, Y = loading) = ' num2str(mean(squeeze(grains_strains_fab_matrix(1, 1, :)))) ', ' num2str(mean(squeeze(grains_strains_fab_matrix(2, 2, :)))) ', ' num2str(mean(squeeze(grains_strains_fab_matrix(3, 3, :))))])
    disp(['Mean strains (surface) (FABLE XYZ, Y = loading) = ' num2str(mean(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) > 300)))) ', ' num2str(mean(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) > 300)))) ', ' num2str(mean(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) > 300))))])
    disp(['standard dev strains (FABLE XYZ, Y = loading) = ' num2str(std(squeeze(grains_strains_fab_matrix(1, 1, :)))) ', ' num2str(std(squeeze(grains_strains_fab_matrix(2, 2, :)))) ', ' num2str(std(squeeze(grains_strains_fab_matrix(3, 3, :))))])
    disp(['standard dev strains (surface) (FABLE XYZ, Y = loading) = ' num2str(std(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) > 300)))) ', ' num2str(std(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) > 300)))) ', ' num2str(std(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) > 300))))])
    %disp(['Mean strains (Kenesei) = ' num2str(mean(grains_strains_ken(:, 1), 1)) ', ' num2str(mean(grains_strains_ken(:, 2), 1)) ', ' num2str(mean(grains_strains_ken(:, 3), 1))])
    disp(' ')
    % Save data
    APS_2016_s1_cycle_1_lattice_strain(end + 1, :) = [mean(squeeze(grains_strains_fab_matrix(1, 1, :))) mean(squeeze(grains_strains_fab_matrix(2, 2, :))) mean(squeeze(grains_strains_fab_matrix(3, 3, :)))];
    APS_2016_s1_cycle_1_lattice_strain_surf(end +1, :) = [mean(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) > 300))) mean(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) > 300))) mean(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) > 300)))];
    APS_2016_s1_cycle_1_lattice_strain_grains{end + 1} = [squeeze(grains_strains_fab_matrix(1, 1, :)) squeeze(grains_strains_fab_matrix(2, 2, :)) squeeze(grains_strains_fab_matrix(3, 3, :))];
    APS_2016_s1_cycle_1_lattice_strain_stdev(end + 1, :) = [std(squeeze(grains_strains_fab_matrix(1, 1, :))) std(squeeze(grains_strains_fab_matrix(2, 2, :))) std(squeeze(grains_strains_fab_matrix(3, 3, :)))];
    APS_2016_s1_cycle_1_lattice_strain_stdev_surf(end +1, :) = [std(squeeze(grains_strains_fab_matrix(1, 1, grains_radial_dist(:) > 300))) std(squeeze(grains_strains_fab_matrix(2, 2, grains_radial_dist(:) > 300))) std(squeeze(grains_strains_fab_matrix(3, 3, grains_radial_dist(:) > 300)))];
    APS_2016_s1_cycle_1_num_of_grains(end +1) = size(grains, 1);
    APS_2016_s1_cycle_1_num_of_grains_surf(end + 1) = sum(grains_radial_dist > 300);
    APS_2016_s1_cycle_1_com{end + 1} = grains_com;
    APS_2016_s1_cycle_1_radius_grains{end + 1} = grains_radius;

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
    APS_2016_s1_cycle_1_orientations{end + 1} = o;
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
        %plotIPDF(o, yvector, 'Property', stiffness_loading'.*squeeze(grains_strains_fab_matrix(2, 2, :)), 'xAxisDirection','east', 'MarkerSize', 4);
        plotIPDF(o, yvector, 'Property', squeeze(grains_strains_fab_matrix(2, 2, :)), 'xAxisDirection','east', 'MarkerSize', 4);
        colorbar; colormap jet; caxis([-3000 5000]);
        print(gcf, ['ipf_strain_step_' num2str(ll)], '-dpng');
        % Visualize COM
        figure;
        scatter3(grains_com(:, 1), grains_com(:, 2), grains_com(:, 3), 121, grains_radial_dist > 400, 'filled');
        campos([0 1800 0]); camtarget([0 200 0]); camup([0 0 1]);
        axis equal;
        axis vis3d;
    end
    close all;
    %clear q r o RMats grains_strains_fab_matrix;
end
break;

%%% FIRST RUN : startup_mtex.m IF IPF COLORBAR NEEDED

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%v
%%% GENERATE IPF COLORMAP USING MTEX
%%% CHECK IN ebsdColorbar.m
% cs  = symmetry('m-3m');
% ss  = symmetry('-1');
% cc  = get_option('antipodal','colorcoding','ipdfHKL');

% [minTheta,maxTheta,minRho,maxRho,v] = getFundamentalRegionPF(cs, 'antipodal');
% h   = S2Grid('PLOT', 'minTheta', minTheta, 'maxTheta', maxTheta,...
%     'minRho', minRho, 'maxRho', maxRho, 'RESTRICT2MINMAX', 'resolution', 1*degree, 'antipodal');
% v   = vector3d(h);
% x   = getx(v); x = x(:);
% y   = gety(v); y = y(:);
% z   = getz(v); z = z(:);
% d   = orientation2color(h,cc,cs,'antipodal');
% save('coloring_scheme.mat', 'x', 'y', 'z', 'd')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%v

% create an EBSD variable containing the data
% DATA
% Sp_ID O[0][0] O[0][1] O[0][2] O[1][0] O[1][1] O[1][2] O[2][0] O[2][1] O[2][2] X Y Z a b c alpha beta gamma Err1 Err2 Err3 MeanRadius Confidence 
% grains	= load('.\ff_hedm\Img66_Layer1_ring4_t300_ring6_t100\Grains.csv'); a0 = 3.014200;    %%%%%%%% PROMISING
% grains	= load('.\ff_hedm\Img66_Layer1_ring4_t100\Grains.csv'); a0 = 3.014200;    %%%%%%%% PROMISING
% grains	= load('.\ff_hedm\Layer1\Grains.csv'); a0 = 3.014200;    %%%%%%%% PROMISING

% grains	= load('.\ff_hedm\Img66_Layer1_ring4_t300_ring6_t100\Grains.csv'); a0 = 3.014200;    %%%%%%%% PROMISING
% grains	= load('.\ff_hedm\Img70_Layer2_ring4_t300_ring6_t100\Grains.csv'); a0 = 3.014200;    %%%%%%%% PROMISING
% grains	= load('.\ff_hedm\Img74_Layer3_ring4_t300_ring6_t100\Grains.csv'); a0 = 3.014200;    %%%%%%%% PROMISING
% grains	= load('.\ff_hedm\Img78_Layer4_ring4_t300_ring6_t100\Grains.csv'); a0 = 3.014200;    %%%%%%%% PROMISING
% grains	= load('.\ff_hedm\Img82_Layer5_ring4_t300_ring6_t100\Grains.csv'); a0 = 3.014200;    %%%%%%%% PROMISING
% grains	= load('.\ff_hedm\Img86_Layer6_ring4_t300_ring6_t100\Grains.csv'); a0 = 3.014200;    %%%%%%%% PROMISING
% grains	= load('.\ff_hedm\Img90_Layer7_ring4_t300_ring6_t100\Grains.csv'); a0 = 3.014200;    %%%%%%%% PROMISING
% grains	= load('.\ff_hedm\Img94_Layer8_ring4_t300_ring6_t100\Grains.csv'); a0 = 3.014200;    %%%%%%%% PROMISING
% grains	= load('.\ff_hedm\Img98_Layer9_ring4_t300_ring6_t100\Grains.csv'); a0 = 3.014200;    %%%%%%%% PROMISING
grains	= load(fullfile('/Volumes/Secondary HD/Work/HEDM_Study_of_Deformation_in_SMA/CHESS_Jun_2015_Steel', 'Grains.csv')); a0 = 3.014200;    %%%%%%%% PROMISING

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FILTER STUFF IF NECESSARY
%%% THRESHOLDING BY COMPLETENESS
Thresh_Completeness = 0.70;
idx_Completeness    = grains(:,24) >= Thresh_Completeness;

%%% THRESHOLDING BY MEAN RADIUS
Thresh_MeanRadius   = 30;
idx_MeanRadius      = grains(:,23) >= Thresh_MeanRadius;

grains  = grains(idx_Completeness & idx_MeanRadius,:);
% grains  = grains(idx_Completeness,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RESRF2APS   = RMatOfQuat(QuatOfESRF2APS);
% RESRF2APS   = eye(3,3);

nGrains     = size(grains, 1);
for i = 1:1:nGrains
    RMats(:,:,i)   =  RESRF2APS*reshape(grains(i,2:10), 3, 3)';
end
qsym    = CubSymmetries; Rsym    = RMatOfQuat(qsym);
quat    = ToFundamentalRegionQ(QuatOfRMat(RMats), qsym);
rod     = RodOfQuat(quat);

xyz = RESRF2APS*[grains(:,11) grains(:,12) grains(:,13)]';
xyz = xyz';
% load('.\coloring_scheme.mat');
% ori     = orientation('quaternion', quat(1,:), quat(2,:), quat(3,:), quat(3,:), cs);
% ebsd    = EBSD(ori, cs, ss);
% rgb     = orientation2color(ori, 'ipdfHSV');
% 
% plot(ebsd, 'colorcoding','ipdfHSV')

%%%% PLOT COM / ONE COLOR
figure, scatter3(xyz(:,1), xyz(:,2), xyz(:,3), grains(:,23)./max(grains(:,23))*100, 'filled', 'b')
grid on
axis square
xlabel('x (um)'); ylabel('y (um)'); zlabel('z (um)')

%%%% PLOT COM / COMPLETENESS AS COLOR
figure, scatter3(xyz(:,1), xyz(:,2), xyz(:,3), 0.01*grains(:,23).*grains(:,23), grains(:,24), 'filled') %% COMPLETENESS
grid on
axis square
xlabel('x (um)'); ylabel('y (um)'); zlabel('z (um)')
return
% %%%% PLOT COM / RGB IN FUNDAMENTAL TRIANGLE AS IPDF
% figure, scatter3(grains(:,11), grains(:,12), grains(:,13), 50, rgb, 'filled') %% COMPLETENESS
% grid on
% axis square

%%%% PLOT ORIENTATIONS / ONE COLOR
figure, PlotFRPerimeter('cubic');
scatter3(rod(1,:), rod(2,:), rod(3,:), 50, 'filled', 'b')
axis square tight off

%%%% PLOT ORIENTATIONS / COMPLETENESS AS COLOR
figure, PlotFRPerimeter('cubic');
scatter3(rod(1,:), rod(2,:), rod(3,:), 64*grains(:,23)*grains(:,23), grains(:,24), 'filled')
axis square tight off
colorbar vert

% %%%% PLOT ORIENTATIONS / RGB IN FUNDAMENTAL TRIANGLE AS IPDF
% figure, PlotFRPerimeter('cubic');
% scatter3(rod(1,:), rod(2,:), rod(3,:), 50, rgb, 'filled')
% axis square tight off

%%%% HISTOGRAM OF GRAIN SIZES NORMALIZED BY MAX GRAIN SIZE
figure, hist(grains(:,23)./max(grains(:,23)), 20)
xlabel('relative grain radius (-)')
ylabel('number of grains (-)')
title(sprintf('Max grain size : %5.0f (micron)', max(grains(:,23))))
% axis([0 1 0 3000])
axis([0 1 0 500])

figure, 
subplot(2,3,1)
hist(grains(:,14))
xlabel('a (Angstrom)')
ylabel('number of grains (-)')
title(sprintf('a0 = %5.4f A', a0))
view([0 90])
% axis([3.00 3.03 0 4000])
axis([3.00 3.03 0 700])
grid on

subplot(2,3,2)
hist(grains(:,15))
xlabel('b (Angstrom)')
ylabel('number of grains (-)')
title(sprintf('a0 = %5.4f A', a0))
view([0 90])
% axis([3.00 3.03 0 4000])
axis([3.00 3.03 0 700])
grid on

subplot(2,3,3)
hist(grains(:,16))
xlabel('c (Angstrom)')
ylabel('number of grains (-)')
title(sprintf('a0 = %5.4f A', a0))
view([0 90])
% axis([3.00 3.03 0 4000])
axis([3.00 3.03 0 700])
grid on

subplot(2,3,4)
hist(grains(:,17))
xlabel('\alpha (degrees)')
ylabel('number of grains (-)')
view([0 90])
% axis([89.7 90.3 0 3000])
axis([89.7 90.3 0 700])
grid on

subplot(2,3,5)
hist(grains(:,18))
xlabel('\beta (degrees)')
ylabel('number of grains (-)')
view([0 90])
% axis([89.7 90.3 0 3000])
axis([89.7 90.3 0 700])
grid on

subplot(2,3,6)
hist(grains(:,19))
xlabel('\gamma (degrees)')
ylabel('number of grains (-)')
view([0 90])
% axis([89.7 90.3 0 3000])
axis([89.7 90.3 0 700])
grid on
