clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MaterialName    = 'Al';                         % FCC Al 
% latticeParms    = 4.050;                        % IN Angstrom
 
% MaterialName    = 'Fe';                         % FCC Fe
% latticeParms    = 3.515;                        % IN Angstrom

MaterialName    = 'Ni';                         % FCC Ni
latticeParms    = 3.520;                        % IN Angstrom

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SampleThickness = 1.0;                          % IN cm
hkls            = load('fcc.hkls');
d_hkls          = PlaneSpacings(latticeParms, 'cubic', hkls');

BeamLineFlux    = load('bm_flux.data');
TakeOffAngle    = 7:1:13;                    % IN deg

for j = 1:1:length(TakeOffAngle)
    wavelength      = 2*d_hkls.*sind(TakeOffAngle(j)/2);
    Energy          = Angstrom2keV(wavelength);         % IN keV
    
    [mu, ~, rho]    = PhotonAttenuation(MaterialName, Energy'./1000);       % cm^2/g
    
    PercentTransmission(:,j)    = exp(-mu*rho*SampleThickness);
    
    Flux    = zeros(length(Energy),1);
    for i = 1:1:length(Energy)
        idx = find(Energy(i) == BeamLineFlux(:,1));
        if isempty(idx)
            idx1    = find(Energy(i) > BeamLineFlux(:,1));
            idx2    = find(Energy(i) < BeamLineFlux(:,1));
            if isempty(idx1) || isempty(idx2)
                Flux(i,1)   = nan;
            else
                idx1    = idx1(end);
                idx2    = idx2(1);
                
                Flux(i,1)   = BeamLineFlux(idx1,2) - (BeamLineFlux(idx1,2) - BeamLineFlux(idx2,2))/(BeamLineFlux(idx1,1) - BeamLineFlux(idx2,1)) * (BeamLineFlux(idx1,1) - Energy(i));
            end
        else
            Flux(i,1)   = BeamLineFlux(idx,2);
        end
    end
    
    TransEnergy(:,j)                    = Energy';
    PhotonTransAtNormalIncidence(:,j)   = PercentTransmission(:,j).*Flux;
end

figure,
set(gcf, 'Position', [1007 33 902 977])

subplot(2,2,1)
semilogy(BeamLineFlux(:,1), BeamLineFlux(:,2), 'k.')
xlabel('Energy (keV)')
ylabel('Number of photons')
title('BM Flux (photons / s / 0.1% BW)')
grid on

subplot(2,2,2)
plot(TransEnergy(:,1), 'b.')
hold on
plot(TransEnergy(:,end), 'r.')
legend(num2str(TakeOffAngle(1)), num2str(TakeOffAngle(end)), 'Location', 'Best')
xlabel('hkl id')
ylabel('Energy (keV)')
title('Diffraction energy - First and the last TOA only')
grid on

subplot(2,2,3)
plot(PercentTransmission(:,1)*100, 'b.')
hold on
plot(PercentTransmission(:,end)*100, 'r.')
axis([1 30 0 100])
legend(num2str(TakeOffAngle(1)), num2str(TakeOffAngle(end)), 'Location', 'Best')
xlabel('hkl id')
ylabel('Percent transmission')
grid on

subplot(2,2,4)
semilogy(PhotonTransAtNormalIncidence(:,1), 'b.')
hold on
semilogy(PhotonTransAtNormalIncidence(:,end), 'r.')
axis([1 30 1e-7 1e14])
legend('7', '13')
xlabel('hkl id')
ylabel('Number of photons transmitted')
grid on