function dydt = phytoplankton_vec(t,y,n,dz,u,D)
P = y;

% set boundary fluxes at 1 + n+1 to 0
J_a(1) = 0;
J_d(1) = 0;
J_a(n+1) = 0;
J_d(n+1) = 0;

i = 2:(n);
% calculate different fluxes
J_a(i) = u*P(i-1);
J_d(i) = -D*(P(i)-P(i-1))/dz;

J = J_a + J_d;

% calculate dPdt (change of phytoplankton in every cell)
dPdt = -(J(2:(n+1))-J(1:n))/dz;
dydt = dPdt';

end