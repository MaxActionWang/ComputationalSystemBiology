%Max Wang, Oct 19, 2021
%course: gene circuit modeling
%accomplished with the help of TA ZYY
%% main function
% The name of the function must be the SAME as the name of the .m script file.
function chain1
global v0; global k1; global k2;

% By changing the following 2 lines, this function simulates transcription activation or inhibition;
%odefunc=@activation;
odefunc=@inhibition;
v0=5;
k1=0.002;
k2=0.0002;

s0=[300,1000,0,0];
tspan = [0,10];

[t,s] = ode45(odefunc, tspan, s0);

figure(1);
plot(t,s(:,1),t,s(:,2),t,s(:,3),t,s(:,4));
xlabel('time')
ylabel('concentration')
legend('X','D','XD','P')
end
%% model function 1: activation
function ds=activation(t,s)
% This function is the model for transcription activation;
% s1~4 stand for the concentration of X, D, XD, P;
global v0; global k1; global k2;
s1=s(1);s2=s(2);s3=s(3);s4=s(4);

ds1 = -k1*s1*s2+k2*s3+v0*s3;
ds2 = -k1*s1*s2+k2*s3;
ds3 = k1*s1*s2-k2*s3-v0*s3;
ds4 = v0*s3;

ds = [ds1;ds2;ds3;ds4];
end
%% model function 2: inhibition
function ds=inhibition(t,s)
%This function is the model of transcription inhibition;
global v0; global k1; global k2;
s1=s(1);s2=s(2);s3=s(3);s4=s(4);
%repressor, X binding region, repressor-DNA, product
ds1 = -k1*s1*s2+k2*s3;
ds2 = -k1*s1*s2+k2*s3;
ds3 = k1*s1*s2-k2*s3;
ds4 = v0*s2/(s2+s3);

ds = [ds1;ds2;ds3;ds4];
end
