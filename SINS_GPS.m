clear all;
clc;
glvs;
glv;
%% ����ԭʼ����
N = 4;                       
Ts = 0.01;                                    
Num=745/Ts;
An = 300/Ts;
GPS_F = 1;
GPS_Fnum = 25;

Data  = load('745s.txt');              
IMU   = Data(:,1:6);                    
Wdata = [IMU(:,1) IMU(:,2) IMU(:,3)];     
Fdata = [IMU(:,4) IMU(:,5) IMU(:,6)];      
GPS   = [Data(:,19) Data(:,20) Data(:,21)]; 
PosR  = [Data(:,13) Data(:,14) Data(:,15)]; 
disp('Step.1:--- ���ݵ������ ---');

%% ��ʼ��׼
AIMU = Data(1:An,:);
[Apsi,Atheta,Agamma] = Align(AIMU);
disp('Step.2:--- ��ʼ��׼���� ---');

%% ��ϵ�����ʼ�� 
Qnb = Trans_att2quat([Apsi Atheta Agamma]*glv.D2R);
Cnb = Trans_att2attm([Apsi Atheta Agamma]*glv.D2R);
[ Vn,Pos,vn_1,vn_2,pos_1,pos_2,posIn_1,posIn_2,Qt,Pk,Xk,Rk,Hk ] = SINS_Init( );

%% ��ϵ���
countNa=1;
countFuse=1;
for count= 1:N:Num-100
    FmIn=Ts*(Fdata(count:count+N-1,:))';   
    WmIn=Ts*(Wdata(count:count+N-1,:))'; 
    %��������
    [Cnb ,Vn, Pos, posIn] = SINS(Cnb,vn_1, vn_2, pos_1, pos_2, posIn_1, posIn_2, WmIn, FmIn, Ts); 
    if rem(countNa,25)==0
        Fn = Cnb*Fdata(countNa*4,:)'; 
        Rm=glv.Re*(1-2*glv.e+3*glv.e*sin(Pos(1))^2);
        Rn=glv.Re*(1+glv.e*sin(Pos(1))^2);  

        Zk=Pos-GPS(countNa*4,:)';
        
        %Kalman�˲�
        [Fhikk_1, Qk] = KalmanPhik(Vn,Cnb,Pos,Fn,Rm,Rn,Qt, GPS_F, 15);
        [Xk, Pk, Kk] = KalmanFilter(Fhikk_1, Qk, Xk, Pk, Hk, Rk, Zk);

        Vn(1:3)=Vn(1:3)-Xk(1:3);
        Pos(1:3)=Pos(1:3)-Xk(7:9);
        Cnb=(eye(3)+antisym_mat(Xk(4:6)))*Cnb;
        
        Na_res(:,countFuse)=[Trans_attm2att(Cnb);Vn;Pos;Xk(10:15);]; 
        Xk(1:9)=0;
        PosT(countFuse,:)=PosR(countNa*4,:)';
        countFuse=countFuse+1;
    end    
    vn_2 = vn_1;                 
    vn_1 = Vn;                  
    pos_2 = pos_1;               
    pos_1 = Pos;                 
    posIn_2 = posIn_1;           
    posIn_1 = posIn;             
    countNa=countNa+1;           
end
disp('Step.3:--- ��ϵ������� ---');

%% ��Ͻ��
figure name '�켣'
plot3(Na_res(7,:)*glv.R2D,Na_res(8,:)*glv.R2D,Na_res(9,:));
hold on ;
plot3(PosT(:,1)'*glv.R2D,PosT(:,2)'*glv.R2D,PosT(:,3));
grid on;
xlabel('γ��/��');ylabel('����/��');zlabel('�߳�/m');

lenRes = [1:length(Na_res)];
figure name 'λ�����'
subplot(311);plot(lenRes,Na_res(7,:)*glv.R2D-PosT(:,1)'*glv.R2D,'r'); 
xlabel('t/s');ylabel('γ�����/��');  title('γ�����');  
grid on; 
subplot(312);plot(lenRes,Na_res(8,:)*glv.R2D-PosT(:,2)'*glv.R2D,'r');
xlabel('t/s');ylabel('�������/��');  title('�������');   
grid on;
subplot(313);plot(lenRes,Na_res(9,:)-PosT(:,3)','r'); 
xlabel('t/s');ylabel('�߳����/m');  title('�߳����'); 
grid on;   