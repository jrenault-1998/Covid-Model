close all
clear all
clc



%model parameters
global alpha C bc ba deltaE deltaIc deltaIa r ...
    deltaSq deltaQ N deltaIp Tf Iclim q0 CT_break s0


alpha = 0.18;             % infection probability
C = 5;                    % contact rate
bc = 0.5;                 % reduction in contacts|symptomatic  (unitless)
ba = 0.75;                % reduction in infectiousness         (unitless)
deltaE = 1/4;             %All "deltaX" terms are (1/days)
deltaIp = 1/3;            %2.4days is right, but needs to be whole number
deltaIc = 1/3.2;          %1/3.2
deltaIa = 1/7;
deltaSq = 1/10;           %How long are people told to isolate for?
deltaQ = 1/10;
ruv = 0.7;                  % (unitless)
N = 5.22e5;                % Pop of NL
q0 = 0.9;                  %CT efficiency
rv = 0.15;                % r for vaccinated individuals

Iclim0 = 0; %Number of symptomatic cases before CT starts
CT_break = 420; %Pop in Ic when CTing breaks down  


%Need to put this a little high
Tf = 180;                    %days of simulation (days)



%%Initial conditions

for z = 1:4

    %4 scenarios for 4 different vaccination status of the pop
    
choice = z;

%α-variant δ-variant
%Efficacy after 1 dose vaccine
%49.2 % (Pfizer) 33.2 % (Pfizer)
%51.4% (AZ) 32.9 % (AZ)
%Efficacy after 2 doses vaccine
%93.4 % (Pfizer) 87.9 % (Pfizer)
%66.1% (AZ) 59.8 % (AZ)

%Assume that we have delta variant, as it is the most prevalent, and Pfizer vaccine.
%0.332 efficiency after 1 dose and 0.879 efficiency after 2 doses ?

if choice == 1
p1 = 0.0;
p2 = 0.75;
s0= (1-p1-p2)*N+p1*N*(1-0.332)+p2*N*(1-0.879);
r = (ruv*s0 + rv*(N-s0))/N;


elseif choice == 2
p1 = 0.0;
p2 = 0.8;
s0= (1-p1-p2)*N+p1*N*(1-0.332)+p2*N*(1-0.879);
r = (ruv*s0 + rv*(N-s0))/N;


elseif choice == 3
p1 = 0.0;
p2 = 0.85;
s0= (1-p1-p2)*N+p1*N*(1-0.332)+p2*N*(1-0.879);
r = (ruv*s0 + rv*(N-s0))/N;


%%% N-s0 %people that cannot be infected or infect.

elseif choice == 4
p1 = 0.0;
p2 = 0.9;
s0= (1-p1-p2)*N+p1*N*(1-0.332)+p2*N*(1-0.879);
r = (ruv*s0 + rv*(N-s0))/N;


end

e0=1;
ip0=0; 
ic0=0;
ia0=0;
Q0=0;
sq0=0;
r0 = 0;
e1=0;
Q1=0; 
Qs=0;
s = [s0; e0; ip0; ic0; ia0; Q0; sq0; r0; e1; Q1;Qs];




%% We plot the number of infected people of the infectious period for different Iclim x efficiency (q0)

%Number of iteractions for variable 1 (Nsteps) and variable 2 (Msteps).


%Step size of the iteractions
dq = 0.15;
dIclim = 2;
%All combinations of all numbers in this vectors will be tested

qin = 0.1;
qvec = qin:dq:1;  

Iclimin = 1;
Iclimvec = Iclimin:dIclim:12;  

Nsteps = size(qvec,2);
Msteps = size(Iclimvec,2);

%Here is the matrix where the final size of the outbreak for each
%combination of the two variables will be saved
Tinf = zeros(Nsteps,Msteps);

for i = 1:Nsteps
for j = 1:Msteps
    
             q0 = qvec(i);
             Iclim = Iclimvec(j);
            
sol = dde23(@CTeq,[1, 2, 3, 4, 5], s,[0 Tf]);

Tnum = sol.y(9,:);
Ic = sol.y(4,:);

%This is just to see how variables change over time at each loop
%figure(3);plot(sol.x,sol.y(9,:),'r');hold on;plot(sol.x,sol.y(4,:),'b');hold on


%Tinf gives the number of infections at the end of the epidemic (or at the
%end of the time interval considered at least) for each combination of the
%two variables of interest

%Set that CT overload if Ic is higher than what CT can handle

if Ic(end) >= CT_break
    Tinf(i,j) = NaN;
else
    Tinf(i,j) = Tnum(end);
end


 end
end

figure(1)
subplot(1,4,z)
%(Put first j variable in x position, then i variable in y position)
h = contourf(Iclimvec,qvec,Tinf);%,'CDataMapping','scaled');
%h = imagesc(Iclimvec,qvec,Tinf,'CDataMapping','scaled');
%h.AlphaData = ones(size(h.CData)); 
%h.AlphaData(isnan(h.CData)) = 0;
colorbar
ylabel('CT efficiency [q_0]')
xlabel('Ic when CT starts [Iclim]')
set(gca,'fontsize',12)
%set(gca,'color',0*[1 1 1]);
set(gca,'YDir','normal')  

end