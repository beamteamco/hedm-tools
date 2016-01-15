clear all
close all
clc

%%% INPUT DATA GENERATED IN STEP0
InputData   = load('example_rsinput.mat');

lambda1 = InputData.lambda1;
lambda2 = InputData.lambda2;
lambda3 = InputData.lambda3;
kappa   = InputData.kappa;

lambda1 = 0;
lambda2 = 0;
lambda3 = 0;
kappa   = 2;

PFNAME_GRID = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_GRID);
GridData    = load(PFNAME_GRID);

PFNAME_MESH = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_MESH);
MeshData    = load(PFNAME_MESH);

plot3(MeshData.x, MeshData.y, MeshData.z, 'bo')
hold on
plot3(GridData.dvc_x, GridData.dvc_y, GridData.dvc_z, 'k^')

% CALC ELLIPSOIDS
[evals,evecs]   = findevs(GridData.dvc_x, GridData.dvc_y, GridData.dvc_z, GridData.ndv);

% evaluates elemental information: shape function derivatives etc at the
% quadrature points
[nqptv, wtq, sfac, ...
    dndxi, dndet, dndze, ...
    nqpts, swt, ssfac, ...
    dnda, dndb]  = shafac( ...
    MeshData.meltyp, MeshData.nnpe, MeshData.nnps);

% Calculate quad point coordinates
[xqpt, yqpt, zqpt]  = qptloc( ...
    MeshData.x, ...
    MeshData.y, ...
    MeshData.z, ...
    MeshData.np,...
    MeshData.numel, ...
    MeshData.nnpe, ...
    sfac);
plot3(xqpt, yqpt, zqpt, 'r.')

% Find the nodal stresses which match the given (or measured) stresses
% while satisfying equilibrium and zero surface traction conditions!
% INITIALIZE VARIABLES
fc  = 1;
kc  = 1;
qc  = 1;

nr_K    = (MeshData.numel*36*MeshData.nnpe*MeshData.nnpe) + ...
    (36*MeshData.numelsyms*MeshData.nnps*MeshData.nnps) + ...
    (36*MeshData.numels*MeshData.nnps*MeshData.nnps);
Ksr = zeros(nr_K, 1);
Ksc = zeros(nr_K, 1);
Ksv = zeros(nr_K, 1);

% MLS METHOD
nr_Q    = MeshData.numel*MeshData.nnpe*36*GridData.ndv;
Qsr     =zeros(nr_Q, 1);
Qsc     =zeros(nr_Q, 1);
Qsv     =zeros(nr_Q, 1);

% loop over the elements to set up matrices
for iele =1:1:MeshData.numel
    [dndx, dndy, dndz, detj]    = sfder( ...
        iele, MeshData.nnpe, nqptv, ...
        sfac, dndxi, dndet, dndze, ...
        MeshData.np, MeshData.x, MeshData.y, MeshData.z);
    
    % Find the weighthing value for the quad points using mls method
    phi = mlsel3d( ...
        xqpt, yqpt, zqpt,...
        GridData.dvc_x, GridData.dvc_y, GridData.dvc_z, ...
        evals, evecs, ...
        kappa, ...
        iele, nqptv, GridData.ndv);
    
    % Residual calculation using mls method
    [se, qe]    = elstif_residual_mls( ...
        MeshData.nnpe, nqptv, GridData.ndv, ...
        wtq, sfac, ...
        dndx, dndy, dndz, ...
        detj, lambda1, phi);
    
    % Assemble K matrix (mls method)
    [Ksr, Ksc, Ksv, kc] = assmbl_ful( ...
        iele, ...
        MeshData.nnpe, MeshData.np, ...
        se, ...
        Ksr, Ksc, Ksv, kc);
    
    % Assemble Q matrix (mls method)
    [Qsr, Qsc, Qsv, qc] = assmbl_ful_q_mls( ...
        iele, MeshData.nnpe, ...
        GridData.ndv, MeshData.np, ...
        qe, Qsr, Qsc, Qsv, qc);
    
    disp(['L2 fit + eq. constraint: ' num2str(iele/MeshData.numel*100) '% is complete!'])
end

% T MATRIX - TRACTION CONDITION
% Calculate the constant matrices
bigNsurf    = bigNsurfmat(MeshData.nnps, nqpts, ssfac);

% For each prescribed free surface
for ieles = 1:1:MeshData.numels
    % Surface normal and jacobian
    [n, rjs]    = surfjac( ...
        ieles, ...
        MeshData.nnps, ...
        nqpts, dnda, dndb, MeshData.nps, ...
        MeshData.x, MeshData.y, MeshData.z);
    
    % Form the surface "stiffness matrix"
    sek = freesurf_residual( ...
        MeshData.nnps, ...
        nqpts, swt, bigNsurf, n, rjs, lambda2);
    
    % Assemble the free surface condition to the corresponding location
    % in the overall stiffness matrix
    [Ksr, Ksc, Ksv, kc]	= assmblsurf_ful( ...
        ieles, ...
        MeshData.nnps, MeshData.nps, ...
        sek, ...
        Ksr, Ksc, Ksv, kc);
    
    disp(['Applying free surface constraint: ' num2str(ieles/MeshData.numels*100) '% is complete!']);
end

% D MATRIX
% symmetry surface condition (in-plane surface tractions are free)
% For each prescribed free surface
for ieles = 1:1:MeshData.numelsyms
    % Surface normal and jacobian
    [n, rjs]    = surfjac( ...
        ieles, ...
        MeshData.nnps, ...
        nqpts, dnda, dndb, MeshData.npsyms, ...
        MeshData.x, MeshData.y, MeshData.z);
    
    % Form the surface "stiffness matrix"
    sek = symsurf_residual( ...
        MeshData.nnps, ...
        nqpts, swt, bigNsurf, n, rjs, lambda3);
    
    % Assemble the free surface condition to the corresponding location
    % in the overall stiffness matrix
    [Ksr, Ksc, Ksv, kc] = assmblsurf_ful( ...
        ieles, ...
        MeshData.nnps, MeshData.npsyms, ...
        sek, ...
        Ksr, Ksc, Ksv, kc);
    
    disp(['Applying symmetry surface constraint: ' num2str(ieles/MeshData.numelsyms*100) '% is complete!'])
end

%%% ASSEMBLE K, Q
K   = sparse(Ksr,Ksc,Ksv);
Q   = sparse(Qsr,Qsc,Qsv);

disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
disp('Generated K and Q matrices')
disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
PFNAME_KQ   = fullfile(InputData.PNAME_SOLUTION, InputData.FNAME_KQ);
disp(['Saving ', InputData.FNAME_KQ])
save(PFNAME_KQ, 'K', 'Q')

pause(1)
disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
disp('Generating Keq')
disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
step6_Generate_matx_Keq

pause(1)
disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
disp('Generating Kfs')
disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
step6_Generate_matx_Kfs

pause(1)
disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
disp('Generating Ksym')
disp('+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
step6_Generate_matx_Ksym