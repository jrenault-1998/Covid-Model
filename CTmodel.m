close all
clear all


%model parameters
global alpha C bc ba k h tau deltaE dektaIr deltaIc deltaIa r epsilon Cv deltasa deltaIv deltaQ deltaQa x1 x2 x3 x4 N deltaIp v

alpha = 1;
C = 1;
bc = 1;
ba = 0.5;
k = 1;
h =1;
tau = 1;
deltaE = 1/4;
dektaIr = 1/2.4;
deltaIc =1/3.2;
deltaIa = 1/7;
deltasa = 1;
deltaIv = 1;
deltaQ = 1;
deltaQa = 1;
deltaIp = 1;
r = 0.7;
epsilon = 1;
Cv = 1;
x1 = 5;
x2 = 6;
x3 = 7;
x4 = 8;
v = 1;
N = 5.22e5;


s0=5.22e5*0.9;
e0=0;
ip=1;
ic=1;
ia=1;
q=0;
qa=0;
sa=0;
sv=5.22e5*0.1;
iv=1;
re=0;
rv=0;

Tf = 120; %days of simulation

options = odeset('RelTol',1e-4,'AbsTol',1e-6);
[T,Y] = ode45(@CTeq, 0:1:Tf, [s0;e0;ip;ic;ia;q;qa;sa;sv;iv;re;rv], options);


figure(1)


subplot(1,2,1)
%title(' ... ')
plot(T,Y(:,2),'k','Linewidth',1)
hold on
plot(T,Y(:,3),'r','Linewidth',1)
hold on
plot(T,Y(:,4),'r','Linewidth',1)
hold on
plot(T,Y(:,5),'r','Linewidth',1)
hold on
plot(T,Y(:,6),'r','Linewidth',1)
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
%legend('I','location','best')
%ylabel('')
xlabel('time [days]')
%axis([ 0 Tf -0.1 55])
set(gca,'fontsize',14)


