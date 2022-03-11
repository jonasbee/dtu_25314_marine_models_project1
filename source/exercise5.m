clear;
close all;
clc;

%% assign parameters and initial values
param.depth = 100;

tspan = 0:1:4000;

param.mu_max = 0.96; % max. growth rate [d^-1]
param.m = 0.24; % mortality [d^-1]
param.u = 1; % sinking speed [m d^-1]
param.D = 2.592; % diffusion coefficient [m^2 d^-1]
param.w = 5; % sinking speed detritus [m d^-1]
param.K_bg = 0.045; % background turbidity [m^-1]
param.k = 5*10^(-10); % light attenuation by phytoplankton [m^2 cell^-1]
param.H_N = 0.0425; % half saturation nutrients [mmol nutrients m^-3]
param.H_I = 20; % half saturation light [µmol photons m^-2 s^-1]
param.alpha = 10^(-9); % nutrient content of phytoplankton [mmol nutrients m^-3]
param.eps = 0.03; % nutrient recycling []
param.gamma = 1.5; % zooplankton grazing []
param.tau = 0.1; % detritus remineralization []
param.N_bot = 75; % nutrients at bottom [mmol nutrients m^-3]

dz_pick = [0.2,0.4,0.8,1.6,3.2,6.4];

% sensitivity analysis parameter
N_bot_pick = [5,20,50,100];
H_N_pick = [0 1 2 3 4 5];

%% vary dz to show independance of grid size
figure
for dz = dz_pick
    param.z = 0.5*dz:dz:(param.depth-(0.5*dz));
    param.n = length(param.z);
    P0 = exp(-(param.z-param.depth/2).^2/5);
    [t,P] = ode45(@phytoplankton_vec,(0:100),P0,[],param.n,dz,param.u,param.D);
    plot(P(end,:), -param.z, 'o-','LineWidth',1)
    drawnow
    hold on
end
title('Varying grid size dz')
legend('dz=0.2','dz=0.4','dz=0.8','dz=1.6','dz=3.2','dz=6.4')
xlabel('phytoplankton concentration (cells/m^3)')
ylabel('depth (m)')

% grid sensitivity shows below 1, grid size has no impact
% assigning 5
param.dz = 1;
param.z = 0.5*param.dz:param.dz:(param.depth-(0.5*param.dz));
param.n = length(param.z);

% inital values
P0 = 10^6*exp(-(param.z-param.depth/2).^2/5);
N0 = param.N_bot*exp(-(param.z-param.depth/1.8).^2/5); 
D0 = 0*exp(-(param.z-param.depth/2).^2/5);
param.PND = [P0 N0 D0];

param.I_0 = 450;

%% perform numerical ode calculation

options = odeset('nonnegative',1:length(param.PND));
[t,y] = ode45(@phytoplankton_li_nu_de_seasonal,tspan,param.PND,options,param);
P = y(:,1:param.n);
N = y(:,(param.n+1):(param.n*2));
D = y(:,((param.n*2)+1):end);

I = seasonal_lightintensity(P(end,:),param,tspan(end));

%% calculate limiting factors
I_lim = I./(param.H_I+I);
N_lim = N(tspan(end),:)./(param.H_N+N(tspan(end),:));


%% local sensitivity analysis (parameters)

% for N_bottom
max_val_N_bot = [];
max_pos_N_bot = [];

zeros(1,length(N_bot_pick));

for N_bot = N_bot_pick
    param.N_bot = N_bot;
    
    [t,y] = ode23(@phytoplankton_li_nu_de_seasonal,tspan,param.PND,options,param);
    P_N_bot = y(:,1:param.n);
    
    m = max_phyto(P_N_bot);
    max_val_N_bot = [max_val_N_bot m(1)];
    max_pos_N_bot = [max_pos_N_bot m(2)];
end

% compartment correction
max_pos_N_bot = max_pos_N_bot.*(param.depth/param.n);

% set H_N back to default
param.N_bot = 75;

% for H_N
max_val_H_N = [];
max_pos_H_N = [];

zeros(1,length(H_N_pick));

for H_N = H_N_pick
    param.H_N = H_N;
    
    [t,y] = ode23(@phytoplankton_li_nu_de_seasonal,tspan,param.PND,options,param);
    P_H_N = y(:,1:param.n);
    
    m = max_phyto(P_H_N);
    max_val_H_N = [max_val_H_N m(1)];
    max_pos_H_N = [max_pos_H_N m(2)];
end

% compartment correction
max_pos_H_N = max_pos_H_N.*(param.depth/param.n);

% set H_N back to default
param.H_N = 0.0425;


%% plotting

% surface plots of P, N, D until converged
figure
subplot(3,1,1)
contourf(t,-param.z,P','EdgeColor', 'none')
colorbar
title('concentration of phytoplankton [cells/m^3]')
ylabel('depth [m]')
xlabel('time [days]')

subplot(3,1,2)
contourf(t,-param.z,N','EdgeColor', 'none')
colorbar
title('concentration of nutrients [mmol nutrient/m^3]')
ylabel('depth [m]')
xlabel('time [days]')

subplot(3,1,3)
contourf(t,-param.z,D','EdgeColor', 'none')
colorbar
title('concentration of detritus [cells/m^3]')
ylabel('depth [m]')
xlabel('time [days]') 

% steady-state solution of model
figure
subplot(1,3,1)
plot(P(end,:),-param.z,'Color','#77AC30','LineWidth',1.5)
title('steady state phytoplankton concentration')
ylabel('depth [m]')
xlabel('phytoplankton concentration [cells/m^3]') 

subplot(1,3,2)
plot(N(end,:),-param.z,'r','LineWidth',1.5)
title('steady state nutrient concentration')
ylabel('depth [m]')
xlabel('nutrient concentration [mmol/m^3]') 

subplot(1,3,3)
plot(D(end,:),-param.z,'-k','LineWidth',1.5)
title('steady state detritus concentration')
ylabel('depth [m]')
xlabel('detritus concentration [cells/m^3]')

% limiting factors
figure
plot(I_lim,-param.z,'Color','#EDB120','LineWidth',1.5)
hold on
plot(N_lim,-param.z,'r','LineWidth',1.5)
plot(P(end,:)./10^8,-param.z,'Color','#77AC30','LineWidth',1.5)
legend('light (I)','nutrients (N)','phytoplankton concentration (P)');
hold off
title('limiting factors')
ylabel('depth [m]')
xlabel('light, nutrients, phytoplankton concentration scaled to 0-1')

% sensitivity of N_bottom
figure
subplot(1,2,1)
plot(N_bot_pick,max_val_N_bot,'LineWidth',1.5)
title('sensitivity of P for N_{bottom}')
xlabel('values for N_{bot}')
ylabel('maximum phytoplankton concentration [cells/m^3]')
subplot(1,2,2)
plot(N_bot_pick,-max_pos_N_bot,'LineWidth',1.5)
title('positions of max. P for N_{bottom} in water column')
xlabel('values for N_{bot}')
ylabel('depth [m]')

% sensitivity of H_N
figure
subplot(1,2,1)
plot(H_N_pick,max_val_H_N,'LineWidth',1.5)
title('sensitivity of P for H_N')
xlabel('values for H_N')
ylabel('maximum phytoplankton concentration [cells/m^3]')
subplot(1,2,2)
plot(H_N_pick,-max_pos_H_N,'LineWidth',1.5)
title('positions of max. P for H_N in water column')
xlabel('values for H_N')
ylabel('depth [m]')

% seasonal behavior of phytoplankton
figure
plot(t,P(:,20)./10^9,'LineWidth',1.5)
hold on
plot(t,N(:,20).*0.5,'LineWidth',1.5)
plot(t,D(:,20)./10^18.5,'LineWidth',1.5)
hold off
title('seasonal succession')
legend('phytoplankton (Px10^{-9})','nutrients (N/2)','detritus (Dx10^{-18.5})')
xlabel('time [days]')
ylabel('scaled concentrations')

figure
plot(t,P(:,20)./10^9,'LineWidth',1.5)
hold on
plot(t,N(:,20).*0.5,'LineWidth',1.5)
plot(t,D(:,20)./10^18.5,'LineWidth',1.5)
xlim([3285 4000])
hold off
title('seasonal succession')
legend('phytoplankton (Px10^{-9})','nutrients (N/2)','detritus (Dx10^{-18.5})')
xlabel('time [days]')
ylabel('scaled concentrations')
