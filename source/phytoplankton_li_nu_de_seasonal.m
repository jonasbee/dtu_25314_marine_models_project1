function dydt = phytoplankton_li_nu_de_seasonal(t,y,param)
%% extract values
P = y(1:param.n)';
N = y(param.n+1:(2*param.n))';
D = y((2*param.n+1):end)';

%% preallocate mu 
mu = zeros(1,param.n);

%% interior fluxes 
i = 2:(param.n);

% phytoplankton
J_a(i) = param.u*(P(i-1));
J_d(i) = (-param.D*(P(i)-P(i-1)))/param.dz;

% diffusion nutrients
J_dN(i) = -param.D.*((N(i)-N(i-1))/param.dz);

% detritus
J_aD(i) = param.w*(D(i-1));
J_dD(i) = (-param.D*(D(i)-D(i-1)))/param.dz;


%% boundary fluxes 
% P
J_a(1) = 0;
J_d(1) = 0;
J_a(param.n+1) = param.u*P(param.n);
J_d(param.n+1) = 0;

J = J_a + J_d;

% N
J_dN(param.n+1) = -param.D*(param.N_bot-N(param.n))/param.dz;
J_dN(1) = 0;

% D
J_aD(1) = 0;
J_dD(1) = 0;
J_aD(param.n+1) = param.w*D(param.n);
J_dD(param.n+1) = 0;

J_D = J_aD + J_dD;

dPdt = -((J(2:(param.n+1))-J(1:param.n))/param.dz);
dNdt = -((J_dN(2:(param.n+1))-J_dN(1:param.n))/param.dz);
dDdt = -((J_D(2:(param.n+1))-J_D(1:param.n))/param.dz);


%% reaction term
I = seasonal_lightintensity(P,param,t);

% growth rate
for i = 1:param.n
    mu(i) = param.mu_max*min((N(i)/(param.H_N+N(i))),(I(i)/(param.H_I+I(i))));
end

% mu(1:param.n) = param.mu_max*min([(N(1:param.n)/(param.H_N+N(1:param.n))),(I(1:param.n)/(param.H_I+I(1:param.n)))]);

%% calculate results
dPdt = dPdt+mu.*P-param.m.*P;
dNdt = dNdt-(param.alpha*mu).*P+param.eps*param.alpha*param.m.*P;
dDdt = dDdt+(param.eps*P+param.gamma*(P.^2)-param.tau*D);

dydt = [dPdt dNdt dDdt]';
end

%% preallocate vector parameters
% J_a = zeros(1,param.n+1);
% J_d = zeros(1,param.n+1);
% J = zeros(1,param.n+1);
% 
% J_dN = zeros(1,param.n+1);
% 
% J_dD = zeros(1,param.n+1);
% J_aD = zeros(1,param.n+1);
% J_D = zeros(1,param.n+1);
% 
% dPdt = zeros(1,param.n);
% dNdt = zeros(1,param.n);
% dDdt = zeros(1,param.n);

% for i = 2:param.n
%     J_a(i) = param.u*(P(i-1));
%     J_d(i) = (-param.D*(P(i)-P(i-1)))/param.dz;
% end
% 
% % diffusion interior fluxes nutrients
% for i = 2:param.n
%     J_dN(i) = -param.D.*((N(i)-N(i-1))/param.dz);
% end
% 
% % interior fluxes for detritus
% for i = 2:param.n
%     J_aD(i) = param.w*(D(i-1));
%     J_dD(i) = (-param.D*(D(i)-D(i-1)))/param.dz;
% end

% for i = 1:param.n
%     dPdt(i) = -((J(i+1)-J(i))/param.dz);
%     dNdt(i) = -((J_dN(i+1)-J_dN(i))/param.dz);
%     dDdt(i) = -((J_D(i+1)-J_D(i))/param.dz);
% end