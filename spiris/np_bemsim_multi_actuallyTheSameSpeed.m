function far_field = np_bemsim_multi(illum,subst,np)
	% % This method didn't turn out to be much faster than np_bemsim, so I'm abandoning it.
% farfieldE = np_bemsim_multi(illum,subst,np)
% 
% Calculate the scattered far fields of an IRIS nanoparticle using
% the boundary element method, implemented in the MNPBEM14 toolbox.
% 
% A structure 
% 
% illum			Structure of properties of each illumination wave 
% 			- Do they all need to have the same wavelength ???
% 
% 	illum{1}.enei 	inverse energy = free space wavelength, nanometers 
% 	illum{2}.pol 	polarization unit vector: e.g. X-polarized: [1 0 0]
% 	illum{3}.dir 	unit vector direction of propagation. 
% 
% subst			Properties of the substrate
% ...
%
% ... etc
% 
% z = 0 is placed at the top of the thin film - eg SiO2-air interface (although a different immersion medium than air can be used)
% 
% far_field is a structure with the scattered fields in response to each excitation

% 1. Set substrate
% ----------------
immersionMaterial 	= epsconst(subst.nList(1)^2);
filmMaterial 		= epsconst(subst.nList(2)^2);
substrateMaterial 	= epsconst(subst.nList(3)^2);
if ischar(np.material) & ( strcmp(np.material,'gold') | strcmp(lower(np.material),'au') )
	particleMaterial = epstable('gold.dat');
else
	particleMaterial = epsconst(np.material^2);
end
materials = {immersionMaterial, filmMaterial, substrateMaterial, particleMaterial};

layerZCoordinates = [0 ,-subst.thickness];
substrate = layerstructure(materials, [1,2,3], layerZCoordinates, layerstructure.options);


% 2. Set particle
% ------------------------
if strcmp(np.type,'sphere')
	% set the np mesh intelligently
	npMeshSize = round(pi*np.size^2/24);
	particle = trisphere( npMeshSize, np.size);
elseif strcmp(np.type,'rod')
	% set the np mesh intelligently
	n1 = round(max([np.size(1), 10]));
	n2 = round(max([np.size(1)/3, 10]));
	n3 = round(np.size(2)/3);
	particle = trirod(np.size(1), np.size(2), [n1,n2,n3]);	particle = rot(particle, 180/pi*acos(np.orientation(3)/norm(np.orientation)), [[0 -1; 1 0]*np.orientation(1:2)'; 0]);
else
	error('np.type must be either ''sphere'' or ''rod''');
end
% np is placed so lowest point is at z = + np.gapToSurface
particle = shift(particle, [ 0, 0, - min( particle.pos( :, 3 ) ) + np.gapToSurface ] );
% Plot the particle
% figure; plot(particle,'nvec', true); axis on;

% 3. Prepare Green's functions
% -----------------
op = bemoptions( 'sim', 'ret',  'interp', 'curv' , 'layer', substrate );
particle = comparticle( materials, { particle }, [ 4, 1 ], 1, op );

tab = tabspace( substrate, particle );
greentab = compgreentablayer( substrate, tab );
op.greentab = greentab;
bem = bemsolver( particle, op );


% 4. Run simulation and measure scattered far fields
% -------------------------------

% pre-set the far field measurement vector
unitSphere = trisphere(2^10,2);
dirVecs = unitSphere.verts;
dirVecs(dirVecs(:,3)<0,:) = []; % We are only interested in back-scattering for IRIS

progressbar('Simulating particle scattering');
for plidx = 1:length(illum)

	% calculate surface charge densities (sig)
	exc = planewave( illum{plidx}.pol, illum{plidx}.dir, op );
	sig = bem \ exc( particle, illum{plidx}.enei );

	spec = spectrum( dirVecs, op );

	f = farfield( spec, sig );

	far_field{plidx}.p = f.p;
	far_field{plidx}.e = f.e;
	far_field{plidx}.enei = f.enei;
	progressbar(plidx/length(illum));
end

end