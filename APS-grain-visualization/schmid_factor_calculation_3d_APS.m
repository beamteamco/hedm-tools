% Schmid factor calculation
schmid_direction = 2; % x = 1, y = 2, z = 3
loading_type = 1;    % -1 = compression, +1 = tension
transformation_system = 2; % 1 = tetragonal, 2 = monoclinic, 3 = Monoclinic (Au)
if(transformation_system == 1)
    % Cubic -> tetragonal
    R  = 0.0;     % Volume change
    es = 0.03;     % Shear
    E = {es*[R-1  0   0; ...
             0   R-1  0; ...
             0    0  R+2 ];
         es*[R-1  0   0; ...
             0   R+2  0; ...
             0    0  R-1 ];
        es*[R+2  0   0; ...
             0   R-1  0; ...
             0    0  R-1 ];
        };
elseif(transformation_system == 2)
    % Cubic -> Monoclinic
    alpha = 1.0243;
    gamma = 0.9563;
    delta = 0.058;
    epsi = -0.0427;
    U = { [gamma  epsi   epsi;      % 1
        epsi   alpha  delta;
        epsi   delta  alpha];
        [gamma -epsi  -epsi;      % 2
        -epsi   alpha  delta;
        -epsi   delta  alpha];
        [gamma -epsi   epsi;      % 3
        -epsi   alpha -delta;
        epsi  -delta  alpha];
        [gamma  epsi  -epsi;      % 4
        epsi   alpha -delta;
        -epsi  -delta  alpha];
        [alpha  epsi   delta;     % 5
        epsi   gamma  epsi;
        delta  epsi  alpha];
        [alpha -epsi   delta;     % 6
        -epsi   gamma -epsi
        delta -epsi  alpha];
        [alpha -epsi  -delta;     % 7
        -epsi   gamma  epsi;
        -delta  epsi  alpha];
        [alpha  epsi  -delta;     % 8
        epsi   gamma -epsi;
        -delta -epsi  alpha];
        [alpha  delta  epsi;      % 9
        delta  alpha  epsi;
        epsi   epsi   gamma];
        [alpha  delta -epsi;      % 10
        delta  alpha -epsi;
        -epsi  -epsi   gamma];
        [alpha -delta  epsi;      % 11
        -delta  alpha -epsi;
        epsi  -epsi   gamma];
        [alpha -delta -epsi;      % 12
        -delta  alpha  epsi;
        -epsi   epsi   gamma];
        };
    E = cell(12,1);
    for i=1:12
        E{i} = 0.5*(transpose(U{i})*U{i} - eye(3));
    end
elseif(transformation_system == 3)
    % Cubic -> monoclinic (James 2014 Nature paper, Au alloy)
    alpha = 1.0591;
    epsi = 0.0073;
    gamma = 0.9363;
    delta = 1.0015;
    U = { [alpha  0   epsi;      % 1
             0   gamma  0;
           epsi   0  delta];
        [alpha  0  -epsi;      % 2
           0  gamma  0;
        -epsi   0  delta];
        [gamma 0   0;      % 3
        0   alpha -epsi;
        0  -epsi  delta];
        [gamma  0  0;      % 4
        0   alpha epsi;
        0  epsi  delta];
        [delta  0   epsi;     % 5
        0   gamma  0;
        epsi  0  alpha];
        [alpha -epsi   0;     % 6
        -epsi   delta 0
        0 0  gamma];
        [alpha epsi  0;     % 7
        epsi   delta  0;
        0  0  gamma];
        [gamma  0  0;     % 8
        0   delta epsi;
        0 epsi  alpha];
        [delta  -epsi  0;      % 9
        -epsi  alpha  0;
        0   0   gamma];
        [gamma  0 0;      % 10
        0  delta -epsi;
        0  -epsi   alpha];
        [delta epsi  0;      % 11
        epsi  alpha 0;
        0  0   gamma];
        [delta 0 -epsi;      % 12
           0  gamma  0;
         -epsi   0   alpha];
        };
    E = cell(12,1);
    for i=1:12
        E{i} = 0.5*(transpose(U{i})*U{i} - eye(3));
    end
end
%
load APS_2016_s1_cycle_1_00_grain_orient_matrices;
%
% In TEXTURE each row is an Euler triplet.
TEXSIZE = size(RMats, 3);
PSCHMID = zeros(TEXSIZE, 1);
TVSCHMID = zeros(TEXSIZE, 1);
TPSCHMID = zeros(TEXSIZE, 1);
% Loop over all Euler triplets.
% Schmid factor for plates
[cb0t, cm0t] = habit_calculation_3d();
for itex=1:TEXSIZE
    %euler = TEXTURE(itex, :)*pi/180.0;
    % Get orientation matrix (crystal to sample!).
    tlg=squeeze(RMats(:, :, itex));
    % Schmid factor for plasticity
    % 12 slip systems. Right now 6 gives [110] systems.
    nslip = 12;
    % Get b, m in crystal coordinate system
    [cb0,cm0] = slipsystem_p(nslip);
    % Transform to sample
    b0 = tlg*cb0;
    m0 = tlg*cm0;
    s0alpha = zeros(3,3,nslip);
    % S = b(x)m
    % Dyadic product is first-column x second-row vector
    % http://en.wikipedia.org/wiki/Dyadic_product
    for isys=1:nslip,
        s0alpha(:, :, isys) = kron(b0(:,isys), m0(:, isys)');
    end
    % Schmid factor for variants
    % Both +ve and -ve Schmid factors are allowed.
    % Hence the abs()
    PSCHMID(itex) = max(abs(s0alpha(schmid_direction,schmid_direction,:)));
    PSCHMIDALL = squeeze(abs(s0alpha(schmid_direction,schmid_direction,:)));
    % Schmid factor for transformation
    % 3 variants
    ntrans = size(E, 1);
    % get b, m in crystal coordinates
    % [cb0,cm0] = slipsystem_t(nslip, 'T');
    % Transform to sample
    Eg = zeros(3,3,ntrans);
    for i=1:ntrans
        Eg(:, :, i) = tlg*E{i}*tlg';
    end
    s0alpha = Eg;
    % -ve Schmid should loose!
    TVSCHMID(itex) = max(loading_type*s0alpha(schmid_direction,schmid_direction,:));
    TVSCHMIDALL = squeeze(loading_type*s0alpha(schmid_direction,schmid_direction,:))';
    nslip = size(cb0t, 2);
    % Transform to sample
    b0=tlg*cb0t;
    m0=tlg*cm0t;
    s0alpha = zeros(3, 3, nslip);
    for isys=1:nslip,
        s0alpha(:, :, isys) = kron(b0(:,isys), m0(:, isys)');
    end
    % -ve Schmid should loose!
    TPSCHMID(itex) = max(loading_type*s0alpha(schmid_direction,schmid_direction,:));
    TPSCHMIDALL = squeeze(loading_type*s0alpha(schmid_direction,schmid_direction,:))';
    %
end

clear Eg TEXSIZE TEX_POLYand b0 m0 cb0 cm0 euler i isys itex m0;
clear ntrans nslip s0alpha tlg;