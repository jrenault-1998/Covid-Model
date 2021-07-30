function dy = CTeq(t,y, Z)

global alpha C bc ba deltaE deltaIp deltaIc deltaIa deltaQ r deltaSq N Iclim q0 CT_break s0

ylag1 = Z(:,1);     %For Z(:,i), these values are t-i days before
ylag2 = Z(:,2);
ylag3 = Z(:,3);
ylag4 = Z(:,4);
ylag5 = Z(:,5);


S  = y(1);
S1 = ylag1(1);
S2 = ylag2(1);      %For "ylagi(j)", this is y(j) at t-i days before
S3 = ylag3(1);
S4 = ylag4(1);
E  = y(2);
E5 = ylag5(2);
Ip = y(3);
Ic = y(4);
Ia = y(5);
Q  = y(6);
Sq = y(7);
R  = y(8);
F = y(9); %Final size of the epidemic
Tq = y(10); %Total size of people in quaratine (infectious)
Tqs = y(11); %Total size of people in quaratine (not infectious)

dy = zeros(11,1);

if (Ic > Iclim)  
    Iclim = 0;
end


steepness = 1000;
shift = pi/2;

function q = qtan(Ic)
      q = q0*(1/pi*(atan(steepness*(Ic-Iclim))+shift));
end


% To see function q
% steepness = 100; 
% shift = pi/2; 
% Iclim = 10;
% CT_break = 420;
% Ic = CT_break:0.1:(CT_break + 100); 
% q = q0*(1/pi*(atan(steepness*(Ic-Iclim))+shift)).*exp(-(1/10)*(Ic - CT_break));
% figure;
% plot(Ic,q)


D = E5*r*deltaE;
q = qtan(Ic);


% Stop running the code if CT capacity is overloaded

if Ic <= CT_break

%S            %Exposed                                    %Contact Traced
dy(1)  = -S*alpha*C*(Ip+bc*Ic+ba*Ia)/N - (D)*q*(1-alpha)*C*(bc*(S+S1)+S2+S3+S4)/N + deltaSq*Sq;

%E            %From S                                     %Contact Traced
dy(2)  = S*alpha*C*(Ip+bc*Ic+ba*Ia)/N - (D)*q*alpha*C*(bc*(S+S1)+S2+S3)/N - deltaE*E;

%Ip       %From E                %Contact Traced        %To Ic
dy(3)  = r*deltaE*E - (D)*r*q*alpha*C*(S4)/N - deltaIp*Ip ;

%Ic                     %To R
dy(4)  = deltaIp*Ip - deltaIc*Ic - q*D;

%Ia        %From E               %Contact Traced                %To R
dy(5)  = (1-r)*deltaE*E - (D)*(1-r)*q*alpha*C*(S4)/N - deltaIa*Ia ;

%Q                 %From E, Ip and Ia                    %To R
dy(6)  = (D)*q*alpha*C*(bc*(S+S1)+S2+S3+S4)/N + q*D - deltaQ*Q;

%Sq                %From S                                    %To S
%%!! Does not include vaccinated individuals put in quarantine!!
dy(7)  = (D)*q*(1-alpha)*C*(bc*(S+S1)+S2+S3+S4)/N - deltaSq*Sq; 

%R       %From Ic      %From Ia    %From Q   %Got Vaccine
dy(8) = deltaIc*Ic + deltaIa*Ia + deltaQ*Q ;

%F %Final size of the epidemic (cumulative number of people getting sick)
dy(9) = deltaE*E;

%Tq %Cumulative number of people in quarantine (infected)
dy(10)  = (D)*q*alpha*C*(bc*(S+S1)+S2+S3+S4)/N + q*D;

%Tqs %Cumulative number of people in quarantine (infected and not)
dy(11)  = (D)*q*(1-alpha)*C*(bc*(S+S1)+S2+S3+S4)/N + ((D)*q*(1-alpha)*C*(bc*(S+S1)+S2+S3+S4)/N)*(N-s0)/S;
%dy(11) = (D)*q*(1-alpha)*C*(bc*(S+S1)+S2+S3+S4)/N * (2-s0/N);



else

dy(1)  = 0;
dy(2)  = 0;
dy(3)  = 0;
dy(4)  = 0;
dy(5)  = 0;
dy(6)  = 0;
dy(7)  = 0; 
dy(8) = 0;
dy(9) = 0;
dy(10)  = 0;

end


end