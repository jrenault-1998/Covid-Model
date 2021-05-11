close all
clear all


%model parameters
global alpha C bc ba h tau deltaE deltaIc deltaIa r epsilon Cv deltasq deltaIv deltaQ deltaQa x1 x2 x3 x4 N deltaIp v

alpha = 0.18;
C = 0.75/alpha;     %beta = alpha*C s.t. beta=0.75
bc = 0.5;           %reduction in contacts|symptomatic?
ba = 0.75;
h =0.8;             %Contact tracing efficacy
tau = 4;            %Estimate, I need to find this ave. value
deltaE = 1/4;
deltaIp = 1/2.4;
deltaIc =1/3.2;
deltaIa = 1/7;
deltasq = 1/10;     %How long are people told to isolate for?
deltaIv = 1/5;      %Average time spent infectious|vaccinated?
deltaQ = 1/10;
deltaQa = 1/10;
r = 0.7;
epsilon = 0.9;
Cv = C + 0.5;       %increase by some constant?
x1 = min(1/deltaIp, tau-1/deltaE);
x2 = min(1/deltaIa, tau);
x3 = x1+r/deltaE;
x4 = x2+(1-r)/deltaE;
v = 0.06/7;         %0.06 of pop. every week
N = 1;
totalpop = 5.2e5;   %Population of Newfoundland


s0=0.9-3/(totalpop);
e0=0.1;
ip=1/(totalpop); %Infectious, Presymptomatic
ic=1/(totalpop);
ia=1/(totalpop);
q=0;
qa=0;
sq=0;
sv=0;
iv=0;
re=0;
rv=0;

Tf = 120; %days of simulation

options = odeset('RelTol',1e-4,'AbsTol',1e-6);
%[T,Y] = ode45(@CTeq, 0:1:Tf, [s0;e0;ip;ic;ia;q;qa;sq;sv;iv;re;rv], options);

options1 = odeset(options,'NonNegative', 1:12);   %Need NonNegative results
[T,Y] = ode45(@CTeq, 0:1:Tf, [s0;e0;ip;ic;ia;q;qa;sq;sv;iv;re;rv], options1);

pop = sum(Y');

figure(1)


subplot(1,2,1)
title('...')
plot(T,Y(:,1),'y','Linewidth',1)
hold on
plot(T,Y(:,2),'k','Linewidth',1)
hold on
plot(T,Y(:,3),'r','Linewidth',1)
hold on
plot(T,Y(:,4),'m','Linewidth',1)
hold on
plot(T,Y(:,5),'c','Linewidth',1)
hold on
plot(T,Y(:,6),'g','Linewidth',1)
%legend('I','location','best')
%ylabel('')
xlabel('time [days]')
%axis([ 0 Tf -0.1 55])
set(gca,'fontsize',14)


subplot(1,2,2)
plot(T,Y(:,7),'r','Linewidth',1)
hold on
plot(T,Y(:,8),'c','Linewidth',1)
hold on
plot(T,Y(:,9),'g','Linewidth',1)
hold on
plot(T,Y(:,10),'b','Linewidth',1)
hold on
plot(T,Y(:,11),'y','Linewidth',1)
%legend('I','location','best')
%ylabel('')
xlabel('time [days]')
%axis([ 0 Tf -0.1 55])
set(gca,'fontsize',14)




