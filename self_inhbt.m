%Max Wang, Oct 19, 2021
%course: gene circuit modeling
% hw4 question 3
%Codes accomplied with the help of TA ZYY
%% main func
function self_inhbt
global alpha;global beta; global k;

%odefunc=@activation;
odefunc=@inhibition;
beta=7;
alpha=0.5;
k=10;

s0=0;
tspan = [0,10];

[t,s] = ode45(odefunc, tspan, s0);

figure(1);
plot(t,s);
xlabel('time')
ylabel('concentration')
legend('X')
end
%% model function
function ds=inhibition(t,s)
global alpha;global beta; global k;
% alpha: degredation rate; beta:maximal production rate; k: inhibition threshold

ds = beta*heaviside(k-s)-alpha*s;
% s= [0:0.01:10];
% plot(s,ds)
end