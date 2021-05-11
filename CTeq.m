function dy = CTeq(t,y)

global alpha C bc ba h tau deltaE deltaIp deltaIc deltaIa deltaQ deltaQa r epsilon Cv x1 x2 x3 x4 v deltasq deltaIv N


S  = y(1);
E  = y(2);
Ip = y(3);
Ic = y(4);
Ia = y(5);
Q  = y(6);
Qa = y(7);
Sq = y(8);
Sv = y(9);
Iv = y(10);
R  = y(11);
Rv = y(12);

dy = zeros(12,1);

dy(1)  = -S*alpha*C*(Ip+bc*Ic+ba*Ia)/N - Ic*h*(1-alpha)*C*tau*S/N + deltasq*Sq - v*N*S/(S+R);
dy(2)  = S*alpha*C*(Ip+bc*Ic+ba*Ia)/N - Ic*h*alpha*C*(1/deltaE)*S/N - deltaE*E;
dy(3)  = r*deltaE*E - Ic*h*alpha*C*x1*S/N - deltaIp*Ip;
dy(4)  = deltaIp*Ip - deltaIc*Ic;
dy(5)  = (1-r)*deltaE*E - Ic*h*alpha*C*x2*S/N - deltaIa*Ia;
dy(6)  = Ic*h*alpha*C*x3*S/N - deltaQ*Q;
dy(7)  = Ic*h*alpha*C*x4*S/N - deltaQa*Qa;
dy(8)  = Ic*h*(1-alpha)*C*tau*S/N - deltasq*Sq; 
dy(9)  = v*N*S/(S+R)-Sv*alpha*Cv*(1-epsilon)*(Ip+bc*Ic+ba*Ia)/N;
dy(10) = Sv*alpha*Cv*(1-epsilon)*(Ip+bc*Ic+ba*Ia)/N - deltaIv*Iv;
dy(11) = deltaIc*Ic+deltaIa*Ia+deltaQ*Q+ deltaQa*Qa-v*N*R/(S+R);
dy(12) = deltaIv*Iv + v*N*R/(S+R);

end