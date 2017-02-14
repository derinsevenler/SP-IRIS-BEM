function [reflected, transmitted] = fresnelFilm(incident, interface)
% Calculates the reflected and transmitted fields after
% a thin film interface is struck by a monochrome plane wave of
% arbitrary polarization.
% The incoming wave should be traveling with a positive z-component.
% The interface has an arbitrary number of layers with real and positive
% refractive indices (no absorption).
% 
% The interface is defined with the first layer starting at z=0.
% 
% 'incident' is a MATLAB structure that describes the incoming wave.
% 
% incident.nu		:	frequency of light in Hz. nu = c/lambda.
% incident.knorm	:	unit 3-vector of propagation direction
% incident.E		:	3-vector - complex magnitude in cartesian coordinates
% The incident wave is defined by
% E_i(r,t) = Re(E*exp(i*(2*pi*nu*t - k[dot]r ) ) )
% 
% 'interface' is a MATLAB structure that describes a structure with N interfaces.
% 
% interface.nList		:	A list of (N+1) indices of refraction, in the order they are encountered by the wave. All must be positive real (no absorption).
% interface.thicknesses	:	A list of (N-1) thicknesses of films, in microns. All must be positive.
% 
% Derin Sevenler - derin@bu.edu, Boston University, 2016.

% Conventions followed are those of Saleh - Introduction to Photonics - Figure 6.2-1
% However, these conventions are internal. The interface simply requires the plane wave
% is traveling in the 

% Assertions:
% 1. Incident wave
% -- real and positive frequency
assert(isreal(incident.nu) & incident.nu > 0, 'Light frequency is not real, positive');
% -- traveling with positive z component - polar angle is between 0 and 90 degrees.
assert(incident.knorm(3)>0, 'Incident wave must have positive z-component');
% -- wave is transverse (plane)
assert(dot(incident.E(:),incident.knorm)<eps, 'Incident wave must be completely transverse (plane)');
% 2. Interface
% -- For N interfaces, there are N+1 indices of refraction and N-1 thicknesses
assert(length(interface.nList) == length(interface.thicknesses)+2, 'Number of refractive indices must be two more than the number of films');
% -- there are at least two refractive indices (no films)
assert(length(interface.nList)>1, 'At least two refractive indices must be provided');
% -- all the indices of refraction are positive and real
assert(isreal(interface.nList) & isempty(find(interface.nList <= 0, 1)), 'Refractive indices must all be positive and real');

% Calculate p-polarized and s-polarized components of the incident field
z_hat = [0 0 1]';
% declination angle theta
theta = zeros(length(interface.nList),1);
theta(1) = acos(incident.knorm(3));
% If the incident wave is normal, let x-direction be s
if theta(1) == 0
	s_hat = [-1 0 0]';
	p_hat = [0 1 0]';
else
	s_hat = cross(z_hat, incident.knorm(:)); s_hat = s_hat/norm(s_hat);
	p_hat = cross(s_hat, incident.knorm(:)); p_hat = p_hat/norm(p_hat);
end
P = dot(incident.E(:), p_hat);
S = dot(incident.E(:), s_hat);

% Calculate the transmission declination angle in each layer
for N = 2:length(interface.nList)
	theta(N) = asin(interface.nList(N-1)/interface.nList(N) * sin(theta(N-1)) );
end

% TODO - Handle internal total reflection (evanescent waves, tunneling?).
assert(isreal(theta), 'There is total internal reflection at some layer!');

% Calculate reflection and transmission coefficients at each interface
nInterfaces = length(interface.nList)-1;

[r_s, r_p, t_s, t_p] = fresnelEquations;
rs = zeros(nInterfaces,1);
rp = zeros(nInterfaces,1);
ts = zeros(nInterfaces,1);
tp = zeros(nInterfaces,1);

for N = 1:nInterfaces
	rs(N) = r_s(interface.nList(N), interface.nList(N+1), theta(N));
	rp(N) = r_p(interface.nList(N), interface.nList(N+1), theta(N));
	ts(N) = t_s(interface.nList(N), interface.nList(N+1), theta(N));
	tp(N) = t_p(interface.nList(N), interface.nList(N+1), theta(N));
end

% Calculate equivalent reflection and transmission coefficients of the whole layered interface
% First interface
Ms = [1, rs(1); rs(1), 1]/ts(1);
Mp = [1, rp(1); rp(1), 1]/tp(1);
% All other interfaces
c = 2.998e14; % speed of light in microns per second
for N = 2:nInterfaces
	% wavelength in this layer in microns
	lambda = c/incident.nu/interface.nList(N);

	% S J Byrnes, http://arxiv.org/pdf/1603.02720v1.pdf
	% Optical path length phase factor
	opl = 2*pi*interface.thicknesses(N-1)*cos(theta(N))/lambda;
	phaseMatrix = [exp(-1i*opl) 0 ; 0 exp(1i*opl) ];
	thisMs = phaseMatrix*[1 rs(N); rs(N) 1]/ts(N);
	thisMp = phaseMatrix*[1 rp(N); rp(N) 1]/tp(N);

	Ms = Ms*thisMs;
	Mp = Mp*thisMp;
end
rs_equiv = Ms(2,1)/Ms(1,1);
ts_equiv = 1/Ms(1,1);
rp_equiv = Mp(2,1)/Mp(1,1);
tp_equiv = 1/Mp(1,1);

% Calculate the reflected and transmitted waves
% Reflected wave - from Snell's Law, x- and y-components of plane wave 
% are conserved, and z-component is flipped.
reflected.knorm = incident.knorm(:).*[1 1 -1]';
% Complex amplitude of the reflected wave
% ... in s- and p-components
s_hat_refl = s_hat;
p_hat_refl = p_hat.*[1 1 -1]';
reflEs = rs_equiv*S*s_hat_refl;
reflEp = rp_equiv*P*p_hat_refl;
% ... in cartesian coordinates
reflected.E = reflEs+reflEp;

% Transmitted wave is rotated by (theta(end) - theta(1)) around s_hat
% https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
thetaRot = theta(end)-theta(1);
transmitted.knorm = incident.knorm(:) * cos(thetaRot) + cross(s_hat, incident.knorm(:))*sin(thetaRot) + s_hat*dot(s_hat, incident.knorm(:))*(1-cos(thetaRot));
% Complex amplitude of transmitted wave
% ... in s- and p-components
s_hat_trans = s_hat;
p_hat_trans = p_hat * cos(thetaRot) + cross(s_hat, p_hat)*sin(thetaRot) + s_hat*dot(s_hat, p_hat)*(1-cos(thetaRot));
transEs = ts_equiv*S*s_hat_trans;
transEp = tp_equiv*P*p_hat_trans;
% ... in cartesian coordinates
transmitted.E = transEs+transEp;

end