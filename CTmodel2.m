close all
clear all
clc

%model parameters
global alpha C bc ba deltaE deltaIc deltaIa totalpop r epsilon1 epsilon2 ...
    Cv deltaSq deltaSv1 deltaIv deltaQ N deltaIp v vmax vstop Tf Iclim tau q0 q02 m %CT_break CT_max 

alpha = 0.18;             %(probability -> unitless)
C = 6;                    %Error for large c and small alpha   (1/day)
bc = 0.5;                 %reduction in contacts|symptomatic?  (unitless)
ba = 0.75;                %reduction in infectiousness         (unitless)
tau = 2;                  %Estimates time from entering Ic to CTing  (days)
deltaE = 1/4;             %All "deltaX" terms are (1/days)
deltaIp = 1/3;            %2.4days is right, but needs to be whole number
deltaIc = 1/3.2;          %1/3.2
deltaIa = 1/7;
deltaSq = 1/10;              %How long are people told to isolate for?
deltaIv = 1/5;               %Average time spent infectious|vaccinated?
deltaQ = 1/10;
deltaSv1 = 1/60;            %Days between first and second dose?
r = 0.7;                    % (unitless)
epsilon1 = 0.6;             % (unitless)
epsilon2 = 0.8;             % (unitless)
Cv = C + 0.5;                %increase by some constant? (1/day)
v = 0.06/7;                  %0.06 of pop. every week  (unitless)
N = 5.2e5;                   %(people)
totalpop = 5.2e5;            %Population of Newfoundland     (people)    
vmax = 462000;               %# of people eligible for vaccine  (unitless)
vstop = 0.5;                 %stop CTing when vstop people are vaccinated (unitless)
q0 = 0.9;
q02 = 0.9;
Iclim = 3;



%Need to put this a little high
Tf = 250;                    %days of simulation (days)

%%I will first leave out these effect
%CT_break = 420;              %Pop in Ic when CTing breaks down  
%CT_max = 450;                %Pop in Ic when CTing bottoms out


%%Initial conditions
s0=5.22e5;
e0=10;
ip0=0; 
ic0=0;
ia0=0;
Q0=0;
sq0=0;
sv1=0;
iv1=0;
iv2=0;
r0=0;
sv2=0;
e1=0;
Q1=0;
s = [s0; e0; ip0; ic0; ia0; Q0; sq0; sv1; iv1; iv2; r0; sv2; e1; Q1];

m = 0;

%% Plots S, Q and Ic in standard 2D lines
sol = dde23(@CTeq,[1, 2, 3, 4, 5], s,[0 Tf]);
figure(1);
plot(sol.x,sol.y(4,:), 'g')
hold on
plot(sol.x,sol.y(3,:)+sol.y(5,:), 'y')
hold on
plot(sol.x,sol.y(13,:), 'm')
hold on
plot(sol.x,(sol.y(6,:)), 'r')
hold on
plot(sol.x,(sol.y(14,:)), 'k')
legend('Ic','Ia + Ip','Cumulative cases','Quarantine','Cumulative Quarantine','Location','Best')
xlabel('time [days]');
ylabel('Population');

figure(2);
plot(sol.x,sol.y(1,:)/N, 'k')
hold on
plot(sol.x,sol.y(8,:)/N, 'g')
hold on
plot(sol.x,sol.y(12,:)/N, 'm')
ylim([0,1]);
legend('S','Sv1','Sv2','Location','Best')
xlabel('time [days]');
ylabel('Population');


