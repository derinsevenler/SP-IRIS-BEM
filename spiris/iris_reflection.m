function reflectedfield = iris_reflection(illum, subst)
% reflectedfield = iris_reflection(illum, subst)
% 
% Calculate the reflected field of the IRIS substrate
% ...
% 

c = 2.998e17; % speed of light in nm/sec
% This function expects illum.dir to be propagating in the -z direction,
% whereas fresnelFilm expects the wave to be traveling in the +z direction.

incident.E = illum.pol.*[1 1 -1];
incident.knorm = illum.dir.*[1 1 -1];
incident.nu = c/illum.enei;

% Input thickness is in nanometers, but fresnelFilm requires microns 
interface.thicknesses = subst.thickness/1e3;
interface.nList = subst.nList;
[reflected, ~] = fresnelFilm(incident, interface);
% Reflect back again
reflectedfield.E = reflected.E.*[1 1 -1]';
reflectedfield.dir = reflected.knorm.*[1 1 -1]';

% %  Use Ronen's reflectivity code to double-check normal incidence
% nVec = subst.nList;
% d = subst.thickness;
% lambda = illum.enei;
% theta = 0;
% rp = 's';
% [r123, t123] = fres2(nVec, d, theta, lambda, rp);

% reflectedfield.E = (illum.pol * r123);
% reflectedfield.dir = illum.dir .* [1 1 -1];

