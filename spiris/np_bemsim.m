function far_field = np_bemsim(illum,subst,np)
% farfieldE = np_bemsim(illum,subst,np)
% 
% Calculate the scattered far fields of an IRIS nanoparticle using
% the boundary element method, implemented in the MNPBEM14 toolbox.
% 
% Excitation: a single plane wave of arbitrary direction and polarization
% 
% illum			Properties of the illumination wave
% 				The illumination is defined as E = Re()
% 
% 	illum.enei 	inverse energy = free space wavelength, nanometers 
% 	illum.pol 	polarization unit vector: e.g. X-polarized: [1 0 0]
% 	illum.dir 	unit vector direction of propagation. 
% 
% subst			Properties of the substrate
% ...
%
% ... etc
% 
% z = 0 is placed at the top of the thin film - eg SiO2-air interface (although a different immersion medium than air can be used)
% (Line 36)

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
	npMeshSize = round(pi*np.size^2/36);
	particle = trisphere( npMeshSize, np.size);
elseif strcmp(np.type,'rod')
	% set the np mesh intelligently
	n1 = round(max([np.size(1), 20]));
	n2 = round(max([np.size(1)/3, 20]));
	n3 = round(max([np.size(2), 15]));
%     n1 = 20;
%     n2 = 20;
%     n3 = 15;
	particle = trirod(np.size(1), np.size(2), [n1,n2,n3]);	particle = rot(particle, 180/pi*acos(np.orientation(3)/norm(np.orientation)), [[0 -1; 1 0]*np.orientation(1:2)'; 0]);
%     figure; plot(particle, 'EdgeColor', 'b' );
%     grid on
%     pause(10);
else
	error('np.type must be either ''sphere'' or ''rod''');
end
% np is placed so lowest point is at z = + np.gapToSurface
particle = shift(particle, [ 0, 0, - min( particle.pos( :, 3 ) ) + np.gapToSurface ] );
% Plot the particle
% figure; plot(particle,'nvec', true); axis on;

% 3. Run simulation
% -----------------
op = bemoptions( 'sim', 'ret',  'interp', 'curv' , 'layer', substrate, 'waitbar', 0);
particle = comparticle( materials, { particle }, [ 4, 1 ], 1, op );

if ~exist( 'greentab', 'var' ) || ~greentab.ismember( substrate, illum.enei, particle )
	tab = tabspace( substrate, particle );
	greentab = compgreentablayer( substrate, tab );
end
op.greentab = greentab;
bem = bemsolver( particle, op );
exc = planewave( illum.pol, illum.dir, op );
sig = bem \ exc( particle, illum.enei );


% 4. Measure scattered far fields
% -------------------------------

unitSphere = trisphere(2^10,2);
dirVecs = unitSphere.verts;
dirVecs(dirVecs(:,3)<0,:) = []; % We are only interested in back-scattering for IRIS

% theta = linspace( 0, 2 * pi, 31 ); % azimuthal
% phi = linspace(0, pi/2, 31); % declination
% thetaPhi = combvec(theta, phi);
% tp = mat2cell(thetaPhi, [2], [repmat(1, 1, size(thetaPhi,2))]);
% sphericalToCartesian = @(x)([sin(x(2)).*cos(x(1)), sin(x(2)).*sin(x(1)), cos(x(2))]);
% dirVcell = cellfun(sphericalToCartesian, tp, 'UniformOutput',false);
% dirVecs= reshape(cell2mat(dirVcell), 3, length(dirVcell))'+1e-6; % I am not sure why I have to add this but it is essential!!

% theta = reshape( linspace( 0, 2 * pi, 50 ), [], 1 );
%  directions for emission
% dirVecs = [ cos( theta ), 0 * theta, sin( theta ) ];

spec = spectrum( dirVecs, op );

f = farfield( spec, sig );
% s = vecnorm( 0.5 * real( cross( f.e, conj( f.h ), 2 ) ) );

far_field.p = f.p;
far_field.e = f.e;
far_field.h = f.h;
far_field.enei = f.enei;


% [ sx, sy ] = pol2cart( theta, s);%30 * s / max( s ) );
% sx = sx*3;
% sy = sy*3;
%  overlay with Poynting vector
% plot( sx, sy, 'w-', 'LineWidth', 1 );

% xx = reshape(dirVecs(:,1).*s, length(theta), length(phi));
% yy = reshape(dirVecs(:,2).*s, length(theta), length(phi));
% zz = reshape(dirVecs(:,3).*s, length(theta), length(phi));
% cc = reshape(s, length(theta), length(phi));
% figure;
% surf(xx,yy,zz, cc);
% axis equal;
