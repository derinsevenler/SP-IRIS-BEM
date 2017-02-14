% Script for testing 'fresnelFilm' function

c = 2.998e14; % speed of light in microns per second
addpath(genpath('/Users/derin/IRIS software/iris-processing/photonics/Refractive Indices'));

lambda0 = linspace(.3,.7,100); % wavelength in microns

interface.nList = [1 0 0]; % initializing first material
interface.thicknesses = [.1]; % microns

theta = 0;%40/57;

incident.E = [cos(theta) 0 sin(theta)]; % x-polarized becomes S-polarized when theta = 0, by convention
incident.knorm = [-sin(theta) 0 cos(theta)]; % normal incidence

for n = 1:length(lambda0)
	
	interface.nList(2) = SiO2RefractiveIndexTemp(lambda0(n), 20);
	interface.nList(3) = SiRefractiveIndexTemp(lambda0(n), 20);

	incident.nu = c/lambda0(n); % frequency of light (Hz)
	[reflected, transmitted] = fresnelFilm(incident, interface);
	refl(n) = norm(reflected.E);
end
figure(2);
plot(lambda0, refl);