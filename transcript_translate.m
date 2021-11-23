%Max Wang, Oct 31, 2021
%course: gene circuit modeling
%hw4, question4
%% main func

global alpha_m;global beta_m; global alpha_p; global beta_p;
alpha_m =0.1; % min^-1
beta_m = 1;% min^-1
alpha_p = 0.05; % min^-1
beta_p = 2; % min^-1

s0=[0,0];
tspan = [0,10];

[t,s] = ode45(@translation, tspan, s0);

figure(1);
plot(t,s);
xlabel('time')
ylabel('concentration')
legend('mRNA','protein')

%% model function
function ds=translation(t,s)
global alpha_m;global beta_m;
global alpha_p;global beta_p;
% alpha: degredation rate; beta:maximal production rate;
s1=s(1); s2=s(2);% mRNA, protein level

ds1 = beta_m-alpha_m*s1;
ds2 = beta_p*s1-alpha_p*s2;

ds = [ds1;ds2];
% s= [0:0.01:10];
% plot(s,ds)
end