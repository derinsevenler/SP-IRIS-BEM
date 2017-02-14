% Script for testing 'fresnelFilm' function

c = 2.998e14; % speed of light in microns per second
addpath(genpath('/Users/derin/IRIS software/iris-processing/photonics/Refractive Indices'));

lambda = .525;
theta = linspace(0, 70*pi/180, 50);
% theta = 45;
% Rays are all in the Y=0 plane - look at knorm below. 
% i.e., phi = 0 in spherical coordinates


interface.nList = [1 SiO2RefractiveIndexTemp(lambda, 20) SiRefractiveIndexTemp(lambda, 20)];
interface.thicknesses = [.1];

incident.nu = c/lambda; % frequency of light (Hz)

for idx = 1:length(theta)
    incident.knorm = [sin(theta(idx)) 0 cos(theta(idx))]; 
	% +z direction, rotated in y=0 plane

    % S-polarized
	incident.E = [0 1 0];
    [reflected, ~] = fresnelFilm(incident, interface);
	refls(idx) = norm(reflected.E);
    refls_ang(idx) = angle(sqrt(sum(reflected.E.^2)));
	
    [r123, ~] = fres2(interface.nList, interface.thicknesses(1), theta(idx), lambda, 's')
    ronenEs = incident.E.*r123; % it's pointed in the same direction
    reflRonenS(idx) = norm(ronenEs);
    reflRonenS_ang(idx) = angle(sqrt(sum(ronenEs.^2)));
    sdif(idx)= sqrt(sum( (reflected.E-transpose(ronenEs)).^2 ) );
    
	% P-polarized
	incident.E = [cos(theta(idx)) 0 -sin(theta(idx))];
    [reflected, ~] = fresnelFilm(incident, interface);
	reflp(idx) = norm(reflected.E);
    reflp_ang(idx) = angle(sqrt(sum(reflected.E.^2)));
    
    [r123, ~] = fres2(interface.nList, interface.thicknesses(1), theta(idx), lambda, 'p')
    ronenEp = incident.E.*r123.*[1 1 -1]; % it has the opposite z-component
    reflRonenP(idx) = norm(ronenEp);
    reflRonenP_ang(idx) = angle(sqrt(sum(ronenEp.^2)));
    pdif(idx)= sqrt(sum( (reflected.E-transpose(ronenEp)).^2 ) );

end
figure;

subplot(3,2,1)
plot(theta*57, refls, theta*57, reflp);
legend('s-polarized','p-polarized');
xlabel('Incident angle');
ylabel('Reflected field amplitude');
title('Reflection Amplitude')

subplot(3,2,2)
plot(theta*57, refls_ang*57, theta*57, reflp_ang*57);
legend('s-polarized','p-polarized');
xlabel('Incident angle');
ylabel('Reflected field phase');
title('Reflection Phase')

subplot(3,2,3)
plot(theta*57, reflRonenS, theta*57, reflRonenP);
legend('s-polarized','p-polarized');
xlabel('Incident angle');
ylabel('Reflected field amplitude');
title('Ronen amplitude')

subplot(3,2,4)
plot(theta*57, reflRonenS_ang*57, theta*57, reflRonenP_ang*57);
legend('s-polarized','p-polarized');
xlabel('Incident angle');
ylabel('Reflected field phase');
title('Ronen Phase')

subplot(3,2,5)
plot(theta*57, abs(sdif), theta*57, abs(pdif));
title('Error')

