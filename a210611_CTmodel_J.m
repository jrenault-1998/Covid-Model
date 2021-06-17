close all
clear all
clc



%model parameters
global alpha C bc ba deltaE deltaIp deltaIc deltaIa deltaQ r deltaSq N m Iclim q0

alpha = 0.18;             % probability of infection given a contact (infection/contact)
C = 5;                    % Contact rate   (contacts/day)
bc = 0.5;                 % reduction in contacts|symptomatic?  (unitless)
ba = 0.75;                % reduction in infectiousness         (unitless)
deltaE = 1/4;             % Time in exposed 
deltaIp = 1/2.4;            % time in presymptomatic (2.4days is right, but needs to be whole number. (why?) )
deltaIc = 1/3.2;          % time as symptomatic (too short?) 1/3.2
deltaIa = 1/7;
deltaSq = 1/10;              %Days of quarantine
deltaQ = 1/10;
r = 0.7;                    % percent of symptomatic 
N = 5.22e5;                      %Tot pop (NL)

Iclim = 5;  %Number of symptomatic cases needed to start contact tracing
q0 = 0.5;  % accuracy of contact tracing

m = 0;

Tf = 50;                    %days of simulation (days)



%% Initial conditions
s0=N;
e0=1;
ip0=0; 
ic0=0;
ia0=0;
Q0=0;
sq0=0;
r0=0;


sol = dde23(@a210611_CTeq_2_J,[1, 2, 3, 4, 5, 6], [s0; e0; ip0; ic0; ia0; Q0; sq0; r0] ,[0 Tf]);

%t = linespace(0,Tf,Tf);
%y = deval(sol,t);

figure(1);
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


        
