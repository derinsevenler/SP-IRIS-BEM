function [r_s, r_p, t_s, t_p] = fresnelEquations(~)
% Fresnel Reflection and Transmission Equations for Jones matrix elements, Saleh 6.2-6, 6.2-7.
% s-polarized is perpendicular to the plane of incidence, in the x-direction in Saleh.
% r is always negative and real.

% Cosine of theta2, from Snell's law
cos2 = @(n1,n2,theta) ( sqrt(1-(n1./n2)^2.*sin(theta).^2) );

r_s = @(n1,n2,theta) ( (n1.*cos(theta) - n2.*cos2(n1,n2,theta))./ (n1.*cos(theta) + n2.*cos2(n1,n2,theta)) );
r_p = @(n1,n2,theta) ( (n1.*cos2(n1,n2,theta) - n2.*cos(theta))./ (n1.*cos2(n1,n2,theta) + n2.*cos(theta)) );

t_s = @(n1,n2,theta) (1 + r_s(n1,n2,theta));
t_p = @(n1,n2,theta) ( (1 + r_p(n1,n2,theta) ).*cos(theta)./cos2(n1,n2,theta) );