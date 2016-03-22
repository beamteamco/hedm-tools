% Plot several grains on IPF
load APS_2016_s1_cycle_1_reconstruction_data;
grain_ids = [527];
load_step = 1;
annotation_property = [527];
%
o = APS_2016_s1_cycle_1_orientations{load_step};
if(exist('annotation_property', 'var'))
    if(size(annotation_property) == size(grain_ids))
        plotIPDF(o(grain_ids), yvector, 'Property', annotation_property, 'xAxisDirection','east', 'MarkerSize', 4);
    end
else
    plotIPDF(o(grain_ids), yvector, 'xAxisDirection','east', 'MarkerSize', 4);
end
