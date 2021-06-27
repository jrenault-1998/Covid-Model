close all
clear all
clc



%model parameters
global alpha C bc ba deltaE deltaIc deltaIa totalpop r epsilon1 epsilon2 Cv deltaSq deltaSv1 deltaIv deltaQ N deltaIp v vmax vstop Tf Iclim tau
global CT_break CT_max

alpha = 0.3;              %(probability -> unitless)
C = 6;                    %Error for large c and small alpha   (1/day)
bc = 0.5;                 %reduction in contacts|symptomatic?  (unitless)
ba = 0.75;                %reduction in infectiousness         (unitless)
Iclim = 5;
tau = 2;                  %Estimates time from entering Ic to CTing  (days)
deltaE = 1/4;             %All "deltaX" terms are (1/days)
deltaIp = 1/3;            %2.4days is right, but needs to be whole number
deltaIc = 1/3.2;          %1/3.2
deltaIa = 1/7;
deltaSq = 1/10;              %How long are people told to isolate for?
deltaIv = 1/5;               %Average time spent infectious|vaccinated?
deltaQ = 1/10;
deltaSv1 = 1/100;            %Days between first and second dose?
r = 0.7;                    % (unitless)
epsilon1 = 0.6;             % (unitless)
epsilon2 = 0.8;             % (unitless)
Cv = C + 0.5;                %increase by some constant? (1/day)
v = 0.06/7;                  %0.06 of pop. every week  (unitless)
N = 5.2e5;                   %(people)
totalpop = 5.2e5;            %Population of Newfoundland     (people)    
vmax = 462000;               %# of people eligible for vaccine  (unitless)
vstop = 0.5;                 %stop CTing when vstop people are vaccinated (unitless)
Tf = 70;                    %days of simulation (days)

CT_break = 420;              %Pop in Ic when CTing breaks down
CT_max = 450;                %Pop in Ic when CTing bottoms out


sol = dde23(@CTeq,[1, 2, 3, 4, 5], @history ,[0 Tf]);


figure(1);
plot(sol.x,sol.y(1,:), 'g')
hold on
plot(sol.x,sol.y(11,:), 'm')
hold on
plot(sol.x,(sol.y(6,:) + sol.y(7,:)), 'r')
hold on
plot(sol.x,(sol.y(8,:) + sol.y(9,:)+ sol.y(10,:)+ sol.y(12,:)), 'y')
legend('S','R','Quarantine','Vaccinated','Location','Best')
xlabel('time [days]');
ylabel('Population');


figure(2);
plot(sol.x,sol.y(3,:),'r:')
hold on
plot(sol.x,sol.y(4,:),'r')
hold on
plot(sol.x,sol.y(5,:),'b')
hold on
plot(sol.x,(sol.y(3,:)+sol.y(4,:)+sol.y(5,:)),'r','linewidth',1.5)
%hold on
%plot(sol.x,sol.y(2,:),'b','linewidth', 1.5)
legend('I_p','I_c','I_a','I_p + I_c + I_a','Location','Best')
xlabel('Time [days]');
ylabel('Active cases');



function s = history(t)
% Constant history function for CTeq.
s0=5.2e5;
e0=0;
ip0=0; %Infectious, Presymptomatic
ic0=1;
ia0=0;
q0=0;
sq0=0;
sv1=0;
iv1=0;
iv2=0;
r0=0;
sv2=0;

s = [s0; e0; ip0; ic0; ia0; q0; sq0; sv1; iv1; iv2; r0; sv2];
end
        

