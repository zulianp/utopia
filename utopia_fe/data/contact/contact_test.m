% This is sphere in contact with infinite halfplane

R1 = 0.5;
R2 = 0;

E1 = 1e2;
E2 = 1e2;

mu1 = 0.3;
mu2 = 0.3;

% d= 0.05;
% d= 0.01;
d = 0.02;
R = R1;

%%%%%%%%%%%%%%%

% E = 1/((1-mu1^2)/E1);

E = 1/((1-mu1^2)/E1 + (1-mu2^2)/E2);


F = (4/3)*E*(R^(1/2))*d^(3/2);

p0= (1/pi)*((6*F*E^2/R^2))^(1/3)

a = (d*R)^(1/2);

r = -0.5:0.01:0.5;

sigma = p0.*(1-(r.^2/a^2)).^(1/2);

close all;
figure;
plot(r, real(sigma));

figure;
plot(r, imag(sigma));