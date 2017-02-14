function filtd = IRISfilt(infields, filtJ)
% Filter some radiated fields with a jones matrix
% 
% unfiltered is a structure: unfiltd.e, unfiltd.p.nlist, etc. same as BEM output.
% 
% filtd is a structure with the same 'fields' as 'unfiltered', but filtered! haha


% -1. Check if it's a scattered or reflected field - they have different notation
if isfield(infields, 'p')
	% It's a scattered field
	unfiltd.nvecs = infields.p.nvec;
	unfiltd.e = infields.e;
else
	% It's a reflected field
	unfiltd.nvecs = infields.dir';
	unfiltd.e = infields.E';
end


% 0. All the fields should already have nvec(:,3)>0, but check just to be sure
collectedRays = find(unfiltd.nvecs(:,3)>0);
Esca = unfiltd.e(collectedRays,:);
nvecs = unfiltd.nvecs(collectedRays,:);

% 1. Find fields in spherical coordinates
theta = atan2(sqrt(unfiltd.nvecs(:,1).^2 + unfiltd.nvecs(:,2).^2 ),unfiltd.nvecs(:,3));
phi = atan2(unfiltd.nvecs(:,2),unfiltd.nvecs(:,1));

nsph_r 		= @(theta, phi) ([sin(theta).*cos(phi),	sin(theta).*sin(phi), 	cos(theta)]);
nsph_theta 	= @(theta, phi) ([cos(theta).*cos(phi),	cos(theta).*sin(phi), 	-sin(theta)]);
nsph_phi 	= @(theta, phi) ([-sin(phi), 			cos(phi), 				zeros(length(phi),1) ]);

Es_theta 	= 	dot(Esca, nsph_theta(theta,phi),2);	% in the azimuthal direction - 'around' the z axis
Es_phi 		= 	dot(Esca, nsph_phi(theta,phi),2); 	% in the declination direction - 'towards' the z axis
Es_r 		= 	dot(Esca, nsph_r(theta,phi),2); 	% in the radial direction - outwards - we expect this field to be zero

% 2. Rotate the fields to cylindrical coordinates.
% We neglect the intensity law, since we would have to un-do it later anyway

Ec_phi 	= Es_phi;
Ec_rho 	= Es_theta;
Ec_z 	= zeros(length(Ec_rho),1);

% 3. Cartesian to cylindrical.

ncyl_phi = @(phi) ([-sin(phi), 			cos(phi), 				zeros(length(phi),1)]); 
ncyl_rho = @(phi) ([cos(phi), 			sin(phi), 				zeros(length(phi),1)]);
ncyl_z = @(phi) ([zeros(length(phi),1),	zeros(length(phi),1), 	ones(length(phi),1)]); % equal for all positions

Ecyl = repmat(Ec_phi,1,3).*ncyl_phi(phi) + repmat(Ec_rho,1,3).*ncyl_rho(phi) + repmat(Ec_z,1,3).*ncyl_z(phi);

% 4. Apply the filter to the fields
filt3 = padarray(filtJ,[1 1], 0,'post');
for n = 1:size(Ecyl,1)
	Ecyl_filtd(n,:) = (filt3*Ecyl(n,:)')';
end

% 5. Back again to cylindrical coordinates!

Efc_phi 	=	dot(Ecyl_filtd, ncyl_phi(phi),2);
Efc_rho		= 	dot(Ecyl_filtd, ncyl_rho(phi),2);
Efc_z		= 	dot(Ecyl_filtd, ncyl_z(phi),2);

% 6. Rotate back to spherical coordinates

Efs_phi 	= Efc_phi;
Efs_theta 	= Efc_rho;
% Efs_r 		= Efs_z; 

% 7. Add these up to cartesian representation
Efs_cart = repmat(Efs_phi,1,3).*nsph_phi(theta,phi) + repmat(Efs_theta,1,3).*nsph_theta(theta,phi);

% package the output correctly
if isfield(infields, 'p')
	filtd.e = Efs_cart;
	filtd.enei = infields.enei;
	filtd.p = struct('nvec', nvecs);
else
	filtd.E = Efs_cart';
	filtd.dir = nvecs';
end