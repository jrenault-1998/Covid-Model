close all
clear all
clc



%model parameters
global alpha C bc ba deltaE deltaIc deltaIa totalpop r epsilon1 epsilon2 ...
    Cv deltaSq deltaSv1 deltaIv deltaQ N deltaIp v vmax vstop Tf Iclim tau q0 q02 %CT_break CT_max 

alpha = 0.18;             %(probability -> unitless)
C = 4;                    %Error for large c and small alpha   (1/day)
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
deltaSv1 = 1/100;            %Days between first and second dose?
r = 0.7;                    % (unitless)
epsilon1 = 0.6;             % (unitless)
epsilon2 = 0.8;             % (unitless)
Cv = C + 0.5;                %increase by some constant? (1/day)
v = 0;%0.06/7;                  %0.06 of pop. every week  (unitless)
N = 5.2e5;                   %(people)
totalpop = 5.2e5;            %Population of Newfoundland     (people)    
vmax = 462000;               %# of people eligible for vaccine  (unitless)
vstop = 0.5;                 %stop CTing when vstop people are vaccinated (unitless)

%Need to put this a little high
Tf = 300;                    %days of simulation (days)

%%I will first leave out these effect
%CT_break = 420;              %Pop in Ic when CTing breaks down  
%CT_max = 450;                %Pop in Ic when CTing bottoms out


%%Initial conditions
s0=5.22e5;
e0=1;
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
s = [s0; e0; ip0; ic0; ia0; Q0; sq0; sv1; iv1; iv2; r0; sv2; e1];

q02 = 0.8; % percentage of symptomatic that contact traces
Iclim = 10;

%% We plot the number of infected people of the infectious period for different q0 and Iclim 

%Number of iteractions for variable 1 (Nsteps) and variable 2 (Msteps).

Nsteps = 10;
Msteps = 10;

%Step size of the iteractions
dq0 = 0.1;
dIclim = 1;

%All combinations of all numbers in this vectors will be tested
q0v = 0:dq0:(Nsteps*dq0);  % accuracy of contact tracing, initial value
q02v = 0:dq0:(Msteps*dq0);

%Here is the matrix where the final size of the outbreak for each
%combination of the two variables will be saved
Tinf = zeros(Nsteps,Msteps);

for i = 1:Nsteps
for j = 1:Msteps
    
    %If you want to change variables, just change the two variables of
    %interest here
             q0 = q0v(i);
            q02 = q02v(j);
   
sol = dde23(@CTeq,[1, 2, 3, 4, 5], s,[0 Tf]);

%Tnum is the vector with the cumulative size of the epidemic at end point
%in time
Tnum = sol.y(13,:);

%Tinf gives the number of infections at the end of the epidemic (or at the
%end of the time interval considered at least) for eah combination of the
%two variables of interest
Tinf(i,j) = Tnum(end);

 end
end

figure(1)
image(q0v,q02v,Tinf)
%colorbar
c = hot(5.22e5);
colormap(c);
%colormap winter
%mycolors = [1 0 0; 1 1 0; 0 0 1];
%colormap(mycolors);
xlabel('q0')
ylabel('q02')
set(gca,'fontsize',14)
set(gca,'YDir','normal') 

%figure(1);
%plot(sol.x,sol.y(1,:), 'g')
% hold on
% plot(sol.x,sol.y(11,:), 'm')
% hold on
% plot(sol.x,(sol.y(6,:) + sol.y(7,:)), 'r')
% hold on
% plot(sol.x,(sol.y(8,:) + sol.y(9,:)+ sol.y(10,:)+ sol.y(12,:)), 'y')
% legend('S','R','Quarantine','Vaccinated','Location','Best')
% xlabel('time [days]');
% ylabel('Population');
% 
% 
% figure(2)
% plot(sol.x,sol.y(3,:),'r:')
% hold on
% plot(sol.x,sol.y(4,:),'r')
% hold on
% plot(sol.x,sol.y(5,:),'b')
% hold on
% plot(sol.x,(sol.y(3,:)+sol.y(4,:)+sol.y(5,:)),'r','linewidth',1.5)
% %hold on
% %plot(sol.x,sol.y(2,:),'b','linewidth', 1.5)
% legend('I_p','I_c','I_a','I_p + I_c + I_a','Location','Best')
% xlabel('Time [days]');
% ylabel('Active cases');


       