function [texture, elements] = makePolyTexture3D_from_MIDAS(number_of_grains, varargin)
% Makes a polycrystalline texture by Voronoi triangulation wrt nuclei
connectivity = importdata('connectivity.inp');
nodes = importdata('nodes.inp');
% Positions of grain nuclei
if(nargin > 2)
    TEX_POLY = varargin{2};
else
    % Orientation of the grains
    if(~exist('TEX_POLY', 'var'))
        load TEX_POLY;
    end
end
if(nargin > 1)
    nuclei = importdata(varargin{1});
    disp(['Reading grain COM from ' varargin{1}])
else
    nuclei = rand(number_of_grains, 3);
end
euler_angles = TEX_POLY(randi(size(TEX_POLY, 1), size(nuclei, 1), 1), :);
elements = zeros(size(connectivity, 1), 3);
connectivity(:, 1) = [];
nodes(:, 1) = [];
for i=1:size(connectivity, 1)
    elements(i, :) = mean(nodes(connectivity(i, :)', :), 1);
end
dist = pdist2(nuclei, elements);
texture = zeros(size(elements, 1), 3);
for i=1:size(connectivity, 1)
    [~, grain_num] = min(dist(:, i), [], 1);
    texture(i, :) = TEX_POLY(grain_num, :);
end
f = fopen('texture.inc', 'w');
fprintf(f, '%4.2f\t %4.2f\t %4.2f\n', texture');
fclose(f);
%
scatter3(elements(:, 1), elements(:, 2), elements(:, 3), 64, texture, 'filled', 'marker', 's');
end