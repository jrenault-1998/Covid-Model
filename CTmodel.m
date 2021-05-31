close all
clear all
clc



%model parameters
global alpha C bc ba q0 deltaE deltaIc deltaIa r epsilon1 epsilon2 Cv deltaSq deltaSv1 deltaIv deltaQ N deltaIp v vmax vstop Tf

alpha = 0.18;
C = 6;                    %Error for large c and small alpha
bc = 0.5;                 %reduction in contacts|symptomatic?
ba = 0.75;
q0 =1;                    %Contact tracing efficacy
tau = 2;                  %Estimates time from entering Ic to CTing
deltaE = 1/4;
deltaIp = 1/3;            %2.4days is right, but needs to be whole number
deltaIc = 1/3.2;             
deltaIa = 1/7;
deltaSq = 1/10;              %How long are people told to isolate for?
deltaIv = 1/5;               %Average time spent infectious|vaccinated?
deltaQ = 1/10;
deltaSv1 = 1/100;            %Days between first and second dose?
r = 0.7;
epsilon1 = 0.6;
epsilon2 = 0.8;
Cv = C + 0.5;                %increase by some constant?
v = 0.06/7;                  %0.06 of pop. every week
N = 1;
totalpop = 5.2e5;            %Population of Newfoundland         
vmax = 462000/totalpop;      %# of people eligible for vaccine
vstop = 0.5;                 %stop CTing when vstop people are vaccinated
Tf = 120;                    %days of simulation


sol = dde23(@CTeq,[1, 2, 3, 4, 5, 6], @history ,[0 Tf]);
figure;
plot(sol.x,sol.y(4,:)*totalpop)
title('Ic');
xlabel('time t');
ylabel('Population');



function s = history(t)
% Constant history function for CTeq.
totalpop = 5.2e5;      %Population of Newfoundland
s0=1-(1/(totalpop));
e0=0;
ip0=1/(totalpop); %Infectious, Presymptomatic
ic0=0;
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
        



% for i = 1:Tf
% %     [TempT,TempY] = ode45(@CTeq, i-1:i, Y(i, :), options);
% %     T(i+1) = TempT(end);
% %     Y(i+1,:) = TempY(end,:);
%     Tempdy = CTeq(i, Y(i, :));
%     
%     Y(i+1, :) = Tempdy' + Y(i, :);
% end

% 
% figure(1)
% 
% 
% subplot(1,2,1)
% title('...')
% plot(T,Y(:,10),'y','Linewidth',1)
% hold on
% plot(T,Y(:,2)+Y(:,3)+Y(:,4)+Y(:,5),'k','Linewidth',1)
% hold on
% plot(T,Y(:,8)+Y(:,9)+Y(:,11),'r','Linewidth',1)
% hold on
% plot(T,Y(:,6)+Y(:,7),'m','Linewidth',1)
% hold on
% plot(T,Y(:,8),'c','Linewidth',1)
% hold on
% plot(T,Y(:,11),'g','Linewidth',1)
% legend('I','location','best')
% ylabel('')
% xlabel('time [days]')
% %axis([ 0 Tf -0.1 55])
% set(gca,'fontsize',14)


% subplot(1,2,2)
% plot(T,Y(:,7),'r','Linewidth',1)
% hold on
% plot(T,Y(:,8),'c','Linewidth',1)
% hold on
% plot(T,Y(:,9),'g','Linewidth',1)
% hold on
% plot(T,Y(:,10),'b','Linewidth',1)
% hold on
% plot(T,Y(:,11),'y','Linewidth',1)
% %legend('I','location','best')
% %ylabel('')
% xlabel('time [days]')
% %axis([ 0 Tf -0.1 55])
% set(gca,'fontsize',14)
% 