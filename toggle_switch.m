% toggle_switch
% Max Wang @Nov 1, 2021
% codes accomplied with the help of TA ZZY
%% main func
R=0;%repressor concentration
I=0;% inducer concentration concentration
s0 = [0,0]; % initial mRNA, protein concentration
tspan=[0,200];
global K; global k_p; global gamma_m; global gamma_p;
gamma_m = 0.19; gamma_p = 2.13*0.01; %min^-1, mRNA and protein degradation
k_m = 2.85; k_p = 3.17; %min^-1, production

D=150; %DNA per cell
K_s = 1000; K_d = 0.05; %repressor
n1=2; n2=2;
K = k_m*D/(1+(R/K_d/(1+(I/K_s)^n2))^n1);

%% part1.1: varied degradation rate of mRNA and protein, with no initial mRNA or protein 
s0 = [0,0]; % initial mRNA, protein concentration
figure()
subplot
subplot(2,1,1)
for gamma_p = [2.13, 4.26, 21.3].*0.01
    gamma_m = 0.19;
    [t,s] = ode45(@half_switch,tspan,s0);
    hold on
    plot(t,s(:,2));
    hold off
end
title('protein level with varied protein degradation rate')
ylabel('protein level')
legend('\gamma_p=2.13*10^{-2} min^{-1}','\gamma_p=4.26*10^{-2} min^{-1}','\gamma_p=2.13*10^{-1} min^{-1}','location','best');

subplot(2,1,2)
for gamma_m = [0.19,0.38,1.9]
    gamma_0 = 0.19;
    [t,s] = ode45(@half_switch,tspan,s0);
    hold on
    plot(t,s(:,1));
    hold off
end
title('mRNA level with varied mRNA degradation rate')
xlabel('time(min)');
ylabel('mRNA level')
legend('\gamma_m=0.19 min^{-1}','\gamma_m=0.38 min^{-1}','\gamma_m=1.9 min^{-1}','location','best');

gamma_m = 0.19; gamma_p = 2.13*0.01;
%% part1.2: varied degradation rate of mRNA and protein, with initial mRNA and protein 
s0 = [2.5*10e3,3.34*10e5]; % initial mRNA, protein concentration
figure()
subplot
subplot(2,1,1)
for gamma_p = [2.13, 4.26, 21.3].*0.01
    gamma_m = 0.19;
    [t,s] = ode45(@half_switch,tspan,s0);
    hold on
    plot(t,s(:,2));
    hold off
end
title('protein level with varied protein degradation rate')
ylabel('protein level')
legend('\gamma_p=2.13*10^{-2} min^{-1}','\gamma_p=4.26*10^{-2} min^{-1}','\gamma_p=2.13*10^{-1} min^{-1}','location','best');

subplot(2,1,2)
for gamma_m = [0.19,0.38,1.9]
    gamma_0 = 0.19;
    [t,s] = ode45(@half_switch,tspan,s0);
    hold on
    plot(t,s(:,1));
    hold off
end
title('mRNA level with varied mRNA degradation rate')
xlabel('time(min)');
ylabel('mRNA level')
legend('\gamma_m=0.19 min^{-1}','\gamma_m=0.38 min^{-1}','\gamma_m=1.9 min^{-1}','location','best');

gamma_m = 0.19; gamma_p = 2.13*0.01;
%% Part 2: equilibrium protein level with varing repressor and protein
[R,I]=meshgrid(0:0.01:5);
K = k_m*D./(1+(R./K_d./(1+(I./K_s).^n2)).^n1);
P_eq = K*k_p/gamma_p/gamma_m+eps;
figure()
mesh(R,I,P_eq)
xlabel('repressor')
ylabel('inducer')
zlabel('protein')
%% Part 3: analytical solution for half switch ODE model

%% model function: half switch
function ds = half_switch(t,s)

global K; global k_p; global gamma_m; global gamma_p;

s1=s(1); s2=s(2);
ds1 = K-gamma_m*s1;
ds2 = k_p*s1-gamma_p*s2;
ds=[ds1;ds2];
end
