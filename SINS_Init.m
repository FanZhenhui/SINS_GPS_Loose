function [ Vn,Pos,vn_1,vn_2,pos_1,pos_2,posIn_1,posIn_2,Qt,Pk,Xk,Rk,Hk ] = SINS_Init( )
glvs; 
Vn = [0 0 0]';                            
Pos = [10*glv.D2R 110*glv.D2R 20]';       

vn_1 = Vn;             
vn_2 = Vn;             
pos_1 = Pos;           
pos_2 = Pos;          
posIn_1 = zeros(3,1);  
posIn_2 = zeros(3,1);  

Qt=zeros(15,15);
Qt(1,1)=(1e-5*glv.G)^2;
Qt(2,2)=(1e-5*glv.G)^2;
Qt(3,3)=(1e-5*glv.G)^2;
Qt(4,4)=(0.001*pi/180/3600)^2;
Qt(5,5)=(0.001*pi/180/3600)^2;
Qt(6,6)=(0.001*pi/180/3600)^2;

Pk=zeros(15,15);       
Pk(1,1)=10^2;          
Pk(2,2)=10^2;         
Pk(3,3)=10^2;          
Pk(4,4)=(5*glv.D2R)^2;           
Pk(5,5)=(5*glv.D2R)^2;           
Pk(6,6)=(5*glv.D2R)^2;  
Pk(7,7)=(1/60*glv.D2R)^2;           
Pk(8,8)=(1/60*glv.D2R)^2;      
Pk(9,9)=(1/60*glv.D2R)^2; 
Pk(10,10)=(10*1e-6*glv.G)^2;           
Pk(11,11)=(10*1e-6*glv.G)^2;      
Pk(12,12)=(10*1e-6*glv.G)^2;       
Pk(13,13)=(0.001*pi/180/3600)^2;         
Pk(14,14)=(0.001*pi/180/3600)^2;           
Pk(15,15)=(0.001*pi/180/3600)^2; 

Xk=zeros(15,1);
Rk=diag([(2.5/glv.Re)^2 (2.5/glv.Re)^2 5^2]);
Hk=[zeros(3,6) eye(3) zeros(3,6)]; 
end