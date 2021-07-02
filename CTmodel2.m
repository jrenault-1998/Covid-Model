close all
clear all
clc

%model parameters
global alpha C bc ba deltaE deltaIc deltaIa totalpop r epsilon1 epsilon2 ...
    Cv deltaSq deltaSv1 deltaIv deltaQ N deltaIp v vmax vstop Tf Iclim tau q0 q02 m %CT_break CT_max 

alpha = 0.18;             %(probability -> unitless)
C = 3;                    %Error for large c and small alpha   (1/day)
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
Tf = 90;                    %days of simulation (days)

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

m = 0;
