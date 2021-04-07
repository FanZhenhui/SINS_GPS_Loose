function [Fhikk_1, Qk] = KalmanPhik(Vn,Cnb,Pos,Fn,Rm,Rn,Qt, Tkf, n)
    glvs;
    F11=[(Vn(2)*tan(Pos(1))-Vn(3))/(Rn+Pos(3)),2*glv.Wie*sin(Pos(1))+(Vn(1)*tan(Pos(1)))/(Rn+Pos(3)),-2*glv.Wie*cos(Pos(1))-Vn(1)/(Rn+Pos(3));
         -2*glv.Wie*sin(Pos(1))-2*(Vn(1)*tan(Pos(1)))/(Rn+Pos(3)),-Vn(3)/(Rm+Pos(3)),-Vn(2)/(Rm+Pos(3));
         2*glv.Wie*cos(Pos(1))+2*Vn(1)/(Rn+Pos(3)),2*Vn(2)/(Rm+Pos(3)),0  ];
    F12=[0,-Fn(3),Fn(2);
         Fn(3),0,-Fn(1);
         -Fn(2),Fn(1),0 ];
    F13=[2*glv.Wie*(Vn(3)*sin(Pos(1))+Vn(2)*cos(Pos(1)))+Vn(1)*Vn(2)*sec(Pos(1))^2/(Rn+Pos(3)) 0 (Vn(1)*Vn(3)-Vn(1)*Vn(2)*tan(Pos(1)))/(Rn+Pos(3))^2;
         -2*Vn(1)*glv.Wie*cos(Pos(1))-Vn(1)^2*sec(Pos(1))^2/(Rn+Pos(3)) 0 Vn(2)*Vn(3)/(Rm+Pos(3))^2+Vn(1)^2*tan(Pos(1))/(Rn+Pos(3))^2;
         -2*Vn(1)*glv.Wie*sin(Pos(1)) 0 -Vn(2)^2/(Rm+Pos(3))^2-Vn(1)^2/(Rn+Pos(3))^2];
    F14=Trans_quat2attm(Cnb);
    F21=[0,-1/(Rm+Pos(3)),0;
         1/(Rn+Pos(3)),0,0;
         tan(Pos(1))/(Rn+Pos(3)),0,0];
    F22=[0,glv.Wie*sin(Pos(1))+Vn(1)*tan(Pos(1))/(Rn+Pos(3)),-glv.Wie*cos(Pos(1))-Vn(1)/(Rn+Pos(3));
         -glv.Wie*sin(Pos(1))-Vn(1)*tan(Pos(1))/(Rn+Pos(3)),0,-Vn(2)/(Rm+Pos(3));
         glv.Wie*cos(Pos(1))+Vn(1)/(Rn+Pos(3)),Vn(2)/(Rm+Pos(3)),0];
    F23=[0 0 Vn(2)/(Rm+Pos(3))^2;
         -glv.Wie*sin(Pos(1)) 0 -Vn(1)/(Rn+Pos(3))^2;
         glv.Wie*cos(Pos(1))+Vn(1)*sec(Pos(1))^2/(Rn+Pos(3)) 0 -Vn(1)*tan(Pos(1))/(Rn+Pos(3))^2];
    F25=-Trans_quat2attm(Cnb);
    F31=[0 1/(Rm+Pos(3)) 0;
         sec(Pos(1))/(Rn+Pos(3)) 0 0;
         0 0 1  ];
    F33=[0 0 -Vn(2)/(Rm+Pos(3))^2;
         Vn(1)*tan(Pos(1))*sec(Pos(1))/(Rn+Pos(3)) 0 -Vn(1)*sec(Pos(1))/(Rn+Pos(3))^2;
         0 0 0  ];
    Ft=[ F11,F12,F13,F14,zeros(3,3);
         F21,F22,F23,zeros(3,3),F25;
         F31,zeros(3,3),F33,zeros(3,6);
         zeros(6,15)]; 
     
    Tkfi = Tkf;     
    facti = 1;      
    Fti = Ft;
    Mi = Qt;
    In = eye(size(Ft,1));
    Fhikk_1 = In + Tkf*Ft;
    Qk = Qt*Tkf;
    for i=2:1:n
        Tkfi = Tkfi*Tkf;        
        facti = facti*i;
        Fti = Fti*Ft;
        Fhikk_1 = Fhikk_1 + Tkfi/facti*Fti;  
        
        FtMi = Ft*Mi;
        Mi = FtMi + FtMi';
        Qk = Qk + Tkfi/facti*Mi;      
    end  
end