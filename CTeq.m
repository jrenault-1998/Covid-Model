function dydt = CTeq(t,y, Z)

global start alpha C bc ba q deltaE deltaIp deltaIc deltaIa deltaQ r epsilon1 epsilon2 Cv  v deltaSq deltaSv1 deltaIv N vmax vstop

ylag1 = Z(:,1);     %For Z(:,i), these values are t-i days before
ylag2 = Z(:,2);
ylag3 = Z(:,3);
ylag4 = Z(:,4);
ylag5 = Z(:,5);
ylag6 = Z(:,6);

S  = y(1);
S1 = ylag1(1);
S2 = ylag2(1);      %For "ylagi(j)", this is y(j) at t-i days before
S3 = ylag3(1);
S4 = ylag4(1);
S5 = ylag5(1);
E  = y(2);
E6 = ylag6(2);
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

<<<<<<< Updated upstream
=======

dy = zeros(12,1);
>>>>>>> Stashed changes


dy = zeros(11,1);

if Sv1 + Iv1 + Iv2 + Sv2 > vmax
    v=0;
<<<<<<< Updated upstream
end


if Sv1 + Iv1 + Iv2 + Sv2 > vstop
    q=0;
=======
else
    v=0.06/7;
>>>>>>> Stashed changes
end

dy(1)  = -S*alpha*C*(Ip+bc*Ic+ba*Ia)/N - (E6*r*deltaE)*q*(1-alpha)*C*(bc*(S+S1)+S2+S3+S4+S5) - v*N*S/(S+R) + deltaSq*Sq;
dy(2)  =  S*alpha*C*(Ip+bc*Ic+ba*Ia)/N - (E6*r*deltaE)*q*alpha*C*(bc*(S+S1)+S2+S3) - deltaE*E;
dy(3)  = r*deltaE*E - (E6*r*deltaE)*r*q*alpha*C*(S4+S5) - deltaIp*Ip;
dy(4)  = deltaIp*Ip - deltaIc*Ic;
<<<<<<< Updated upstream
dy(5)  = (1-r)*deltaE*E - (E6*r*deltaE)*(1-r)*q*alpha*C*(S4+S5) - deltaIa*Ia;
dy(6)  = (E6*r*deltaE)*q*alpha*C*(bc*(S+S1)+S2+S3+S4+S5) - deltaQ*Q;
dy(7)  = (E6*r*deltaE)*q*(1-alpha)*C*(bc*(S+S1)+S2+S3+S4+S5) - deltaSq*Sq; 
dy(8)  = v*N*S/(S+R)-Sv1*alpha*Cv*(1-epsilon1)*(Ip+bc*Ic+ba*Ia)/N - deltaSv1*Sv1;
dy(9)  = Sv1*alpha*Cv*(1-epsilon1)*(Ip+bc*Ic+ba*Ia)/N - deltaIv*Iv1;
dy(10) = Sv2*alpha*Cv*(1-epsilon2)*(Ip+bc*Ic+ba*Ia)/N - deltaIv*Iv1;
dy(11) = deltaIc*Ic + deltaIa*Ia + deltaQ*Q - v*N*R/(S+R);
dy(12) = v*N*R/(S+R) + deltaIv*Iv1 + deltaSv1*Sv1 + deltaIv*Iv2 - Sv2*alpha*Cv*(1-epsilon2)*(Ip+bc*Ic+ba*Ia)/N;
=======
dy(5)  = (1-r)*deltaE*E - Ic*h*alpha*C*x2*S/N - deltaIa*Ia;
dy(6)  = Ic*h*alpha*C*x3*S/N - deltaQ*Q;
dy(7)  = Ic*h*alpha*C*x4*S/N - deltaQa*Qa;
dy(8)  = Ic*h*(1-alpha)*C*tau*S/N - deltasq*Sq; 
dy(9)  = v*N*S/(S+R)-Sv*alpha*Cv*(1-epsilon)*(Ip+bc*Ic+ba*Ia)/N;
dy(10) = Sv*alpha*Cv*(1-epsilon)*(Ip+bc*Ic+ba*Ia)/N - deltaIv*Iv;
dy(11) = deltaIc*Ic+deltaIa*Ia+deltaQ*Q+ deltaQa*Qa-v*N*R/(S+R);
dy(12) = deltaIv*Iv + v*N*R/(S+R);
>>>>>>> Stashed changes
