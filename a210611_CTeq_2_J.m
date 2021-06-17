function dy = a210611_CTeq_2_J(t,y, Z)

global alpha C bc ba deltaE deltaIp deltaIc deltaIa deltaQ r deltaSq N m Iclim q0 q 

ylag1 = Z(:,1);     %For Z(:,i), these values are t-i days before
ylag2 = Z(:,2);
ylag3 = Z(:,3);
ylag4 = Z(:,4);
ylag5 = Z(:,5);
%ylag6 = Z(:,6);
%ylag7 = Z(:,7);
%ylag8 = Z(:,8);


S  = y(1);
S1 = ylag1(1);
S2 = ylag2(1);      %For "ylagi(j)", this is y(j) at t-i days before
S3 = ylag3(1);
S4 = ylag4(1);
E  = y(2);
E5 = ylag5(2);
Ip = y(3);
Ic = y(4);
%Ic1 = ylag1(4);
Ia = y(5);
Q  = y(6);
Sq = y(7);
R  = y(8);



dy = zeros(8,1);


if (Ic > Iclim)  
    q = q0;
else
    q = 0;
end

if (Ic > Iclim)  
    Iclim = 0.1;
end




%S            %Exposed                                    %Contact Traced
dy(1)  = -S*alpha*C*(Ip+bc*Ic+ba*Ia)/N - (E5*r*deltaE)*q*(1-alpha)*C*(bc*(S+S1)+S2+S3+S4)/N + deltaSq*Sq;

%E            %From S                                     %Contact Traced
dy(2)  = S*alpha*C*(Ip+bc*Ic+ba*Ia)/N - (E5*r*deltaE)*q*alpha*C*(bc*(S+S1)+S2+S3)/N - deltaE*E;

%Ip       %From E                %Contact Traced        %To Ic
dy(3)  = r*deltaE*E - (E5*r*deltaE)*r*q*alpha*C*(S4)/N - deltaIp*Ip + r*m;

%Ic                     %To R
dy(4)  = deltaIp*Ip - deltaIc*Ic - E5*r*deltaE*q;

%Ia        %From E               %Contact Traced                %To R
dy(5)  = (1-r)*deltaE*E - (E5*r*deltaE)*(1-r)*q*alpha*C*(S4)/N - deltaIa*Ia + (1-r)*m;

%Q                 %From E, Ip and Ia                    %To R
dy(6)  = (E5*r*deltaE)*q*alpha*C*(bc*(S+S1)+S2+S3+S4)/N +Q*Ic- deltaQ*Q + q*E5*r*deltaE;

%Sq                %From S                                    %To S
dy(7)  = (E5*r*deltaE)*q*(1-alpha)*C*(bc*(S+S1)+S2+S3+S4)/N - deltaSq*Sq; 

%R       %From Ic      %From Ia    %From Q   %Got Vaccine


dy(8) = deltaIc*Ic + deltaIa*Ia + deltaQ*Q;




end