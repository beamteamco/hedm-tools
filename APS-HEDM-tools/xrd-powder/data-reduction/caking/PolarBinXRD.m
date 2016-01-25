function polImg = PolarBinXRD(mesh, instr, cakeParms, img)
% Polar integration of image data

% clear all
% close all
% clc
% load('PolarBinXRD_input.mat')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
img   = double(img);

L   = instr.detectorsize/instr.pixelsize;

% !!! THESE ARE IN THE CARTESIAN FRAME !!!
x0  = cakeParms.origin(1);   % in pixels
y0  = cakeParms.origin(2);   % in pixels

startRho    = cakeParms.sector(3);
stopRho     = cakeParms.sector(4);

numAzi  = cakeParms.bins(1);             %%% NUMBER OF AZIMUTHAL POINTS
numRho  = cakeParms.bins(2);             %%% NUMBER OF RADIAL POINTS PER AZIMHUTH
numEta  = cakeParms.bins(3);             %%% NUMBER OF ETA POINTS PER AZIMUTH
dAzi    = 360/numAzi;
dEta    = dAzi/numEta;

R   = cakeParms.sector(3):(cakeParms.sector(4)-cakeParms.sector(3))/numRho:cakeParms.sector(4);
R   = repmat(R, numEta + 1, 1)';

Rlist   = mean((R(1:end - 1,:) + R(2:end,:))./2,2);

polImg.azimuth   = cakeParms.azim;
polImg.radius    = zeros(numRho, numAzi);
polImg.intensity = zeros(numRho, numAzi);

for ii = 1:1:numAzi
    fprintf('Processing sector %d of %d\n', ii, numAzi);
    tic;
    
    azi_ini = polImg.azimuth(ii) - dAzi/2;
    azi_fin = polImg.azimuth(ii) + dAzi/2;
    TH  = azi_ini:dEta:azi_fin;
    TH  = repmat(TH, numRho + 1, 1);
    
    [x, y]  = pol2cart(deg2rad(TH),R);
    x       = x0 + x;
    y       = y0 + y;
    
    V   = zeros(numRho + 1, numEta + 1);
    for i = 1:1:(numRho + 1)
        for j = 1:1:(numEta + 1)
            xy  = [x(i,j); y(i,j)];
            
            V(i,j) = DataCoordinates(xy, L, mesh, img);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SOME BENCHMARK TESTS
            % gridNum     = (L - 1)*(fix(xy(2)) - 1) + fix(xy(1));
            % elemNums    = [gridNum*2-1 gridNum*2];
            % 
            % P1  = mesh.crd(:, mesh.con(:,elemNums(1)));
            % P2  = mesh.crd(:, mesh.con(:,elemNums(2)));
            % IN1 = inpolygon(xy(1), xy(2), P1(1,:), P1(2,:));
            % IN2 = inpolygon(xy(1), xy(2), P2(1,:), P2(2,:));
            % if IN1
            %     A1  = [P1(1,1)-P1(1,3) P1(1,2)-P1(1,3);P1(2,1)-P1(2,3) P1(2,2)-P1(2,3)];
            %     B1  = [xy(1)-P1(1,3); xy(2)-P1(2,3)];
            %     Y1  = A1\B1;
            %     
            %     fele    = elemNums(1);
            %     fcrd    = [Y1(1) Y1(2) 1-Y1(1)-Y1(2)];
            % elseif IN2
            %     A2  = [P2(1,1)-P2(1,3) P2(1,2)-P2(1,3);P2(2,1)-P2(2,3) P2(2,2)-P2(2,3)];
            %     B2  = [xy(1)-P2(1,3); xy(2)-P2(2,3)];
            %     Y2  = A2\B2;
            %     
            %     fele    = elemNums(2);
            %     fcrd    = [Y2(1) Y2(2) 1-Y2(1)-Y2(2)];
            % else
            %     fprintf('UH-OH')
            % end
            % fcon    = img(mesh.con(:, fele));
            % V(i,j)  = dot(fcon, fcrd', 1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    
    if cakeParms.fastint
        V       = mean(V,2);
        V       = (V(1:end-1) + V(2:end))/2;
        Ilist   = V;
    else
        Ilist   = BuildMeshPolarXRD(R, V, mesh.qrule);
    end
    
    polImg.radius(:,ii)    = Rlist';
    polImg.intensity(:,ii) = Ilist';
    
    dtime   = toc;
    fprintf('Processing time for sector %d is %1.4f\n', ii, dtime);
    % return
end
% return
polImg.radius       = polImg.radius';
polImg.intensity    = polImg.intensity';