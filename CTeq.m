function dy = CTeq(t,y, Z)

global alpha C bc ba deltaE deltaIp deltaIc deltaIa deltaQ totalpop r epsilon1 epsilon2 Cv  v deltaSq deltaSv1 deltaIv N vmax vstop Iclim
global CT_break CT_max


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
Sv1 = y(8);
Iv1 = y(9);
Iv2 = y(10);
R  = y(11);
Sv2 = y(12);



dy = zeros(11,1);



if Sv1 + Iv1 + Iv2 + Sv2 > vmax             %Vaccination stops at %population eligible
    v=0;
end


steepness = 5;
shift = 1/2;

function Iclim0 = Iclimit(Ic)
  if Ic > 20/totalpop
      Iclim0 = 0.1;
      
  else 
      
      Iclim0 = Iclim;
  end
end

Iclim = Iclimit(Ic);

function q = qtan(Ic)
  if Ic*totalpop < CT_break
      q = (1/pi)*(atan(steepness*(Ic*totalpop-Iclim)))+shift;
      
  else
      q = ((1/pi)*(atan(steepness*(Ic*totalpop-Iclim)))+shift)*exp(-(Ic*totalpop - CT_break));
      
  end
end



q = qtan(Ic);

disp(q);

% 
% if Sv1 + Iv1 + Iv2 + Sv2 > vstop            %Contact tracing stops at %population vaccinated
%     q = 0;
% else
%     
% end


%S            %Exposed                                    %Contact Traced
dy(1)  = -S*alpha*C*(Ip+bc*Ic+ba*Ia)/N - (E5*r*deltaE)*q*(1-alpha)*C*(bc*(S+S1)+S2+S3+S4) - v*N*S/(S+R) + deltaSq*Sq;

%E            %From S                                     %Contact Traced
dy(2)  = S*alpha*C*(Ip+bc*Ic+ba*Ia)/N - (E5*r*deltaE)*q*alpha*C*(bc*(S+S1)+S2+S3) - deltaE*E;

%Ip       %From E                %Contact Traced        %To Ic
dy(3)  = r*deltaE*E - (E5*r*deltaE)*r*q*alpha*C*(S4) - deltaIp*Ip;

%Ic                     %To R
dy(4)  = deltaIp*Ip - deltaIc*Ic;

%Ia        %From E               %Contact Traced                %To R
dy(5)  = (1-r)*deltaE*E - (E5*r*deltaE)*(1-r)*q*alpha*C*(S4) - deltaIa*Ia;

%Q                 %From E, Ip and Ia                    %To R
dy(6)  = (E5*r*deltaE)*q*alpha*C*(bc*(S+S1)+S2+S3+S4) - deltaQ*Q;

%Sq                %From S                                    %To S
dy(7)  = (E5*r*deltaE)*q*(1-alpha)*C*(bc*(S+S1)+S2+S3+S4) - deltaSq*Sq; 

%Sv1     %From S                 %To Iv1                            %To Sv2
dy(8)  = v*N*S/(S+R) - Sv1*alpha*Cv*(1-epsilon1)*(Ip+bc*Ic+ba*Ia)/N - deltaSv1*Sv1;

%Iv1                %From Sv1                            %To Sv2
dy(9)  = Sv1*alpha*Cv*(1-epsilon1)*(Ip+bc*Ic+ba*Ia)/N - deltaIv*Iv1;

%Iv2                %From Sv2                            %To Sv2
dy(10) = Sv2*alpha*Cv*(1-epsilon2)*(Ip+bc*Ic+ba*Ia)/N - deltaIv*Iv2;

%R       %From Ic      %From Ia    %From Q   %Got Vaccine
dy(11) = deltaIc*Ic + deltaIa*Ia + deltaQ*Q - v*N*R/(S+R);

%Sv2      %------Vaccine & Recovered--------        %2nd Dose     %To Iv2
dy(12) = v*N*R/(S+R) + deltaIv*Iv1 + deltaIv*Iv2 + deltaSv1*Sv1 - Sv2*alpha*Cv*(1-epsilon2)*(Ip+bc*Ic+ba*Ia)/N;




end