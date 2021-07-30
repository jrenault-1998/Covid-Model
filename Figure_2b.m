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
r = 0.7;                  % (unitless)
N = 5.22e5;                % Pop of NL
q0 = 0.9;                  %CT efficiency

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
%0.332 efficiency after 1 dose and 0.879 efficacy after 2 doses ?

if choice == 1
p1 = 0.0;
p2 = 0.75;
s0= (1-p1-p2)*N+p1*N*(1-0.332)+p2*N*(1-0.879);

elseif choice == 2
p1 = 0.0;
p2 = 0.8;
s0= (1-p1-p2)*N+p1*N*(1-0.332)+p2*N*(1-0.879);

elseif choice == 3
p1 = 0.0;
p2 = 0.85;
s0= (1-p1-p2)*N+p1*N*(1-0.332)+p2*N*(1-0.879);

%%% N-s0 %people that cannot be infected or infect.

elseif choice == 4
p1 = 0.0;
p2 = 0.9;
s0= (1-p1-p2)*N+p1*N*(1-0.332)+p2*N*(1-0.879);

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


%% We plot the number of infected people of the infectious period for different CT efficiency(q) x importation rate(m)

%Number of iteractions for variable 1 (Nsteps) and variable 2 (Msteps).


%Step size of the iteractions
dq0 = 0.15;
dm = 1;
%All combinations of all numbers in this vectors will be tested

q0in = 0.1;
q0v = q0in:dq0:1;  % accuracy of contact tracing, initial value


mmin = 1;
mmax = 7;
mvec = mmin:dm:mmax;  % infectiousness of COVID, initial value

Nsteps = size(q0v,2);
Msteps = size(mvec,2);


%Here is the matrix where the final size of the outbreak for each
%combination of the two variables will be saved
Tinf = zeros(Nsteps,Msteps);

for i = 1:Nsteps
for j = 1:Msteps
    
             q0 = q0v(i);
             m = mvec(j);
            
    Iclim = Iclim0;


if m == 1

sol1 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[0 Tf]);

%set fix timepoints
tint1 = 0:1:Tf;
sol1f = deval(sol1,tint1);

%Define Tnum1 and Ic1
Tnum1 = sol1f(9,:);
Ic1 = sol1f(4,:);


Ic = Ic1;
Tnum = Tnum1;


elseif m == 2

Iclim = Iclim0;
sol2 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf/2)  Tf]);

%set fix timepoints
tint2 = round(Tf/2):1:Tf;
sol2f = deval(sol2,tint2);


Tnum2p = sol2f(9,:);
Ic2p = sol2f(4,:);

Tnum2 = [zeros(1,round(Tf/2))  Tnum2p];
Ic2 = [zeros(1,round(Tf/2))  Ic2p];

Ic = Ic1+Ic2;
Tnum = Tnum1+Tnum2;

elseif m == 3

    Iclim = Iclim0;
    sol3 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf/3) Tf]);
    Iclim = Iclim0;
    sol32 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf*2/3) Tf]);

%set fix timepoints
tint3 = round(Tf/3):1:Tf;
sol3f = deval(sol3,tint3);

tint32 = round(Tf*2/3):1:Tf;
sol32f = deval(sol32,tint32);


Tnum3p = sol3f(9,:);
Ic3p = sol3f(4,:);

Tnum32p = sol32f(9,:);
Ic32p = sol32f(4,:);


Tnum3 = [zeros(1,round(Tf/3))  Tnum3p];
Ic3 = [zeros(1,round(Tf/3))  Ic3p];

Tnum32 = [zeros(1,round(Tf*2/3))  Tnum32p];
Ic32 = [zeros(1,round(Tf*2/3))  Ic32p];


Ic = Ic1+Ic32+Ic3;
Tnum = Tnum1+Tnum32+Tnum3;


elseif m == 4

sol4 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf/4) Tf]);
Iclim = Iclim0;
sol42 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf*2/4) Tf]);
Iclim = Iclim0;
sol43 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf*3/4) Tf]);

%set fix timepoints
tint4 = round(Tf/4):1:Tf;
sol4f = deval(sol4,tint4);

tint42 = round(Tf*2/4):1:Tf;
sol42f = deval(sol42,tint42);
tint43 = round(Tf*3/4):1:Tf;
sol43f = deval(sol43,tint43);

Tnum4p = sol4f(9,:);
Ic4p = sol4f(4,:);

Tnum42p = sol42f(9,:);
Ic42p = sol42f(4,:);
Tnum43p = sol43f(9,:);
Ic43p = sol43f(4,:);

Tnum4 = [zeros(1,round(Tf/4))  Tnum4p];
Ic4 = [zeros(1,round(Tf/4))  Ic4p];

Tnum42 = [zeros(1,round(Tf*2/4))  Tnum42p];
Ic42 = [zeros(1,round(Tf*2/4))  Ic42p];
Tnum43 = [zeros(1,round(Tf*3/4))  Tnum43p];
Ic43 = [zeros(1,round(Tf*3/4))  Ic43p];

Ic = Ic1+Ic42+Ic43+Ic4;
Tnum = Tnum1+Tnum42+Tnum43+Tnum4;

elseif m == 5

    Iclim = Iclim0;
sol5 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf/5) Tf]);
Iclim = Iclim0;
sol52 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf*2/5) Tf]);
Iclim = Iclim0;
sol53 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf*3/5) Tf]);
Iclim = Iclim0;
sol54 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf*4/5) Tf]);

%set fix timepoints
tint5 = round(Tf/5):1:Tf;
sol5f = deval(sol5,tint5);

tint52 = round(Tf*2/5):1:Tf;
sol52f = deval(sol52,tint52);
tint53 = round(Tf*3/5):1:Tf;
sol53f = deval(sol53,tint53);
tint54 = round(Tf*4/5):1:Tf;
sol54f = deval(sol54,tint54);


Tnum5p = sol5f(9,:);
Ic5p = sol5f(4,:);

Tnum52p = sol52f(9,:);
Ic52p = sol52f(4,:);
Tnum53p = sol53f(9,:);
Ic53p = sol53f(4,:);
Tnum54p = sol54f(9,:);
Ic54p = sol54f(4,:);

Tnum5 = [zeros(1,round(Tf/5))  Tnum5p];
Ic5 = [zeros(1,round(Tf/5))  Ic5p];

Tnum52 = [zeros(1,round(Tf*2/5))  Tnum52p];
Ic52 = [zeros(1,round(Tf*2/5))  Ic52p];
Tnum53 = [zeros(1,round(Tf*3/5))  Tnum53p];
Ic53 = [zeros(1,round(Tf*3/5))  Ic53p];
Tnum54 = [zeros(1,round(Tf*4/5))  Tnum54p];
Ic54 = [zeros(1,round(Tf*4/5))  Ic54p];

Ic = Ic1+Ic52+Ic53+Ic54+Ic5;
Tnum = Tnum1+Tnum52+Tnum53+Tnum54+Tnum5;


elseif m == 6

    Iclim = Iclim0;
sol6 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf/6) Tf]);
Iclim = Iclim0;
sol62 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf*2/6) Tf]);
Iclim = Iclim0;
sol63 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf*3/6) Tf]);
Iclim = Iclim0;
sol64 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf*4/6) Tf]);
Iclim = Iclim0;
sol65 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf*5/6) Tf]);


%set fix timepoints
tint6 = round(Tf/6):1:Tf;
sol6f = deval(sol6,tint6);

tint62 = round(Tf*2/6):1:Tf;
sol62f = deval(sol62,tint62);
tint63 = round(Tf*3/6):1:Tf;
sol63f = deval(sol63,tint63);
tint64 = round(Tf*4/6):1:Tf;
sol64f = deval(sol64,tint64);
tint65 = round(Tf*5/6):1:Tf;
sol65f = deval(sol65,tint65);

Tnum6p = sol6f(9,:);
Ic6p = sol6f(6,:);

Tnum62p = sol62f(9,:);
Ic62p = sol62f(6,:);
Tnum63p = sol63f(9,:);
Ic63p = sol63f(6,:);
Tnum64p = sol64f(9,:);
Ic64p = sol64f(6,:);
Tnum65p = sol65f(9,:);
Ic65p = sol65f(6,:);


Tnum6 = [zeros(1,round(Tf/6))  Tnum6p];
Ic6 = [zeros(1,round(Tf/6))  Ic6p];


Tnum62 = [zeros(1,round(Tf*2/6))  Tnum62p];
Ic62 = [zeros(1,round(Tf*2/6))  Ic62p];
Tnum63 = [zeros(1,round(Tf*3/6))  Tnum63p];
Ic63 = [zeros(1,round(Tf*3/6))  Ic63p];
Tnum64 = [zeros(1,round(Tf*4/6))  Tnum64p];
Ic64 = [zeros(1,round(Tf*4/6))  Ic64p];
Tnum65 = [zeros(1,round(Tf*5/6))  Tnum65p];
Ic65 = [zeros(1,round(Tf*5/6))  Ic65p];

Ic = Ic1+Ic62+Ic63+Ic64+Ic65+Ic6;
Tnum = Tnum1+Tnum62+Tnum63+Tnum64+Tnum65+Tnum6;


elseif m == 7
    

Iclim = Iclim0;
sol7 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf/7) Tf]);
Iclim = Iclim0;
sol72 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf*2/7) Tf]);
Iclim = Iclim0;
sol73 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf*3/7) Tf]);
Iclim = Iclim0;
sol74 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf*4/7) Tf]);
Iclim = Iclim0;
sol75 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf*5/7) Tf]);
Iclim = Iclim0;
sol76 = dde23(@CTeq,[1, 2, 3, 4, 5], s,[round(Tf*6/7) Tf]);

%set fix timepoints
tint7 = round(Tf/7):1:Tf;
sol7f = deval(sol7,tint7);

tint72 = round(Tf*2/7):1:Tf;
sol72f = deval(sol72,tint72);
tint73 = round(Tf*3/7):1:Tf;
sol73f = deval(sol73,tint73);
tint74 = round(Tf*4/7):1:Tf;
sol74f = deval(sol74,tint74);
tint75 = round(Tf*5/7):1:Tf;
sol75f = deval(sol75,tint75);
tint76 = round(Tf*6/7):1:Tf;
sol76f = deval(sol76,tint76);


Tnum7p = sol7f(9,:);
Ic7p = sol7f(6,:);

Tnum72p = sol72f(9,:);
Ic72p = sol72f(6,:);
Tnum73p = sol73f(9,:);
Ic73p = sol73f(6,:);
Tnum74p = sol74f(9,:);
Ic74p = sol74f(6,:);
Tnum75p = sol75f(9,:);
Ic75p = sol75f(6,:);
Tnum76p = sol76f(9,:);
Ic76p = sol76f(6,:);


Tnum7 = [zeros(1,round(Tf/7))  Tnum7p];
Ic7 = [zeros(1,round(Tf/7))  Ic7p];

Tnum72 = [zeros(1,round(Tf*2/7))  Tnum72p];
Ic72 = [zeros(1,round(Tf*2/7))  Ic72p];
Tnum73 = [zeros(1,round(Tf*3/7))  Tnum73p];
Ic73 = [zeros(1,round(Tf*3/7))  Ic73p];
Tnum74 = [zeros(1,round(Tf*4/7))  Tnum74p];
Ic74 = [zeros(1,round(Tf*4/7))  Ic74p];
Tnum75 = [zeros(1,round(Tf*5/7))  Tnum75p];
Ic75 = [zeros(1,round(Tf*5/7))  Ic75p];
Tnum76 = [zeros(1,round(Tf*6/7))  Tnum76p];
Ic76 = [zeros(1,round(Tf*6/7))  Ic76p];


Ic = Ic1+Ic72+Ic73+Ic74+Ic75+Ic76+Ic7;
Tnum = Tnum1+Tnum72+Tnum73+Tnum74+Tnum75+Tnum76+Tnum7;


end

%Tinf gives the number of infections at the end of the epidemic (or at the
%end of the time interval considered at least) for each combination of the
%two variables of interest

%Set that CT overload if Ic is higher than what CT can handle

if Ic(end) >= CT_break
    Tinf(i,j) = NaN;
else
    Tinf(i,j) = Tnum(end);
end

%Plot to see the timeseries at each interaction
%figure(3);plot(0:1:Tf,Ic,'r');hold on;plot([0:1:Tf],Tnum,'b');hold on


end




end


figure(1)
%(Put first j variable in x position, then i variable in y position)
subplot(1,4,z)
%h = imagesc(mvec,q0v,Tinf,'CDataMapping','scaled');
h = contourf(mvec,q0v,Tinf); %,'CDataMapping','scaled');
%h.AlphaData = ones(size(h.CData)); 
%h.AlphaData(isnan(h.CData)) = 0;
colorbar
ylabel('CT efficiency [q_0]')
xlabel('Importation rate [m]')
set(gca,'fontsize',12)
set(gca,'YDir','normal')  

end