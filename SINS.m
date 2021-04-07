function [cnb,vn,pos,dRnk]=SINS(cnb_1,vn_1,vn_2,pos_1,pos_2,dRnk_1,dRnk_2,wm,vm,ts)
glvs;
Tm=4*ts;
%% ×ËÌ¬
Lk = 1.5*pos_1(1)-0.5*pos_2(1);
hk = 1.5*pos_1(3)-0.5*pos_2(3); 
Wniek = [0;glv.Wie*cos(Lk);glv.Wie*sin(Lk)]; 
wnen = [-vn_1(2)/(glv.Re+hk);vn_1(1)/(glv.Re+hk);vn_1(1)*tan(Lk)/(glv.Re+hk)];
Fck = [0,-1/(glv.Re+hk),0;1/(glv.Re+hk),0,0;tan(Lk)/(glv.Re+hk),0,0];
Cge = [-sin(pos_1(2)),cos(pos_1(2)),0;-sin(pos_1(1))*cos(pos_1(2)),-sin(pos_1(1))*sin(pos_1(2)),cos(pos_1(1));cos(pos_1(1))*cos(pos_1(2)),cos(pos_1(1))*sin(pos_1(2)),sin(pos_1(1))];
Sk = Wniek*Tm+Fck*(2*dRnk_1-dRnk_2);
att_nb_1=Trans_attm2att(cnb_1);
qnb_1   =Trans_att2quat(att_nb_1);
phi=wm(:,1)+wm(:,2)+wm(:,3)+wm(:,4)+736/945*(cross(wm(:,1),wm(:,2))+cross(wm(:,3),wm(:,4)))+334/945*(cross(wm(:,1),wm(:,3))+cross(wm(:,2),wm(:,4)))+526/945*cross(wm(:,1),wm(:,4))+654/945*cross(wm(:,2),wm(:,3));
psi_k=Wniek*Tm+Fck*(2*dRnk_1-dRnk_2);
Cn=eye(3)+sin(norm(psi_k))/norm(psi_k)*antisym_mat(psi_k)+(1-cos(norm(psi_k)))/(norm(psi_k))^2*antisym_mat(psi_k)*antisym_mat(psi_k);
Cb=eye(3)+sin(norm(phi))/norm(phi)*antisym_mat(phi)+(1-cos(norm(phi)))/(norm(phi))^2*antisym_mat(phi)*antisym_mat(phi);
Cnb=Cn'*cnb_1*Cb;
cnb = Cnb;

%% ËÙ¶È
vnk=1.5*vn_1-0.5*vn_2;
gk=glv.G;
gk=[0;0;-gk];
dvngcork=gk*Tm-cross((Fck*vnk+2*Wniek),(2*dRnk_1)-dRnk_2);
dvbk = sum(vm,2)+0.5*cross(sum(wm,2),sum(vm,2))+cross((54/105*wm(:,1)+92/105*wm(:,2)+214/105*wm(:,3)),vm(:,4))+cross((54/105*vm(:,1)+92/105*vm(:,2)+214/105*vm(:,3)),wm(:,4));
dvnsfk=quat_mulVec(qnb_1,dvbk)-0.5*antisym_mat(Sk)*Trans_quat2attm(qnb_1)*sum(vm,2);
vn=vn_1+dvnsfk+dvngcork;

%% Î»ÖÃ
dthet1=wm(:,1)+wm(:,2);
dthet2=wm(:,3)+wm(:,4);
dv1=vm(:,1)+vm(:,2);
dv2=vm(:,3)+vm(:,4);
A=1/Tm*(3*dv1-dv2);
B=4/Tm^2*(dv2-dv1);
Sdvk=Tm^2*A/2+Tm^3*B/6;
dRrotk=Tm*(cross(dthet1,(5/18*dv1+1/6*dv2))+cross(dthet2,(1/6*dv1+1/18*dv2)));
dRsculk=Tm*(cross(dthet1,(11/90*dv1+1/10*dv2))+cross(dthet2,(1/90*dv2-7/30*dv1)));
dRbsfkk_1=Sdvk+dRrotk+dRsculk;
dRnsfk=-1/3*antisym_mat(Sk)*quat_mulVec(qnb_1,dvbk)*Tm+quat_mulVec(qnb_1,dRbsfkk_1);
dRnk=(vn_1+0.5*dvngcork)*Tm+dRnsfk;
zetak=Fck*dRnk;
cne=eye(3)-antisym_mat(zetak);
cpos=cne*Cge;
cpos=cpos*(cpos'*cpos)^(-1/2);
Latitude=asin(cpos(3,3))*180/pi;
longmain=atan(cpos(3,2)/cpos(3,1))*180/pi;
if longmain<0 && cpos(3,1)<0
    Longtitude=longmain+180;
end
if longmain>0 && cpos(3,1)<0
    Longtitude=longmain-180;
end
if cpos(3,1)>0
    Longtitude=longmain;
end
h=pos_1(3)+dRnk(3);
pos=[Latitude*pi/180;Longtitude*pi/180;h];