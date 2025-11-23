%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab demo code for the paper
% A novel Gaussian-Student¡¯s t-Skew mixture distribution based Kalman filter
% H Zou, S Wu, Q Xue, X Sun, M Li
% Signal Processing


% If you use our code in your publication, please cite the above paper.
% Demo code written by Han Zou.
% Email: 932501425@qq.com
% Looking forward to your feedback.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

profile on;
clear all;
close all;
clc;
randn('state',sum(100*clock));
format long;

%%%%Model parameters
nxp=1;
nx=4;
nz=2;
T=1;
q=0.1;
r=2;

F=[eye(2) T*eye(2);zeros(2) eye(2)];
H=[eye(2) zeros(2)];
Q0=[T^3/3*eye(2) T^2/2*eye(2);T^2/2*eye(2) T*eye(2)]*q;
R0=r*eye(2);
ts=400;

%%%%Parameter selections
U1=10;      %%%%%%Parameter of state outlier
U2=100;      %%%%%%Parameter of measurement outlier  USELESS
N=100;        %%%%%%Maximum number of iterations


for expt = 1:nxp
    
    fprintf('MC number = %d\n',expt);
    
    %%%%Initial values
    x=[0;0;10;10];
    P=diag([100 100 100 100]);
    Skk=utchol(P);

    xf=x+Skk*randn(nx,1);
    Pf=P;
    beta1_bar=1.5*ones(nx,1);
    sigama1=1;
    x_proposed=xf;
    P_proposed=Pf;
    v_stgt_1=5;
    v_stgt_2=5;
    c10=0.8;
    c20=0.1;
    c30=0.1;
    g10=0.8;
    g20=0.1;
    g30=0.1;
    f0=5;
    u0=5;
    beta2_bar=1.5*ones(nz,1);
    sigama2=1;
    tao_p=5;
    tao_r=5;
    xA=x;
    xpro=xf;
    for t=1:ts
        
        %%%%%%%Simulate non-stationary noise
        %%%%%%%Gaussian noise
        if t<=100
            p1=0.8;
            p2=0.8;
        end
        %%%%%%%Slightly heavy-tailed noise
        if (t>100)&&(t<=200)
            p1=0.6;
            p2=0.6;
        end
        %%%%%%%Moderately heavy-tailed noise
        if (t>200)&&(t<=300)
            p1=0.4;
            p2=0.4;
        end
        %%%%%%%Gaussian noise
        if t>300
            p1=0.2;
            p2=0.2;
        end
        
        %%%%Simulate true state and measurement
        test1=rand;
        test2=rand;
        
        if test1<=p1
            Q=Q0;
            wk=utchol(Q)*randn(nx,1);%%¸ßË¹ÔëÉù
        else
            Q=100*Q0;
            wk=utchol(Q)*randn(nx,1)+skew_noise(diag([1.0;1.0;1.0;1.0]),Q0,5)+utchol(Q0)*randn(nx,1);
        end
        
        R=10*R0;
        
        if test2<=p2 
            vk=utchol(R0)*randn(nz,1);
        else
            vk=utchol(R)*randn(nz,1);
        end
        
        %%%%True noise covariance matrices
        D_Q=p1*Q0+(1-p1)*U1*Q0;
        D_R=R0;
        
        %%%%Simulate true state and measurement
        x=F*x+wk;
        z=H*x+vk;
            
        %%%%Filtering

        [x_proposed,P_proposed]=proposed_algorithm(x_proposed,P_proposed,F,H,z,Q0,R0,...
            f0,u0,c10,c20,c30,g10,g20,g30,tao_p,tao_r,v_stgt_1,v_stgt_2,beta1_bar,sigama1,beta2_bar,sigama2,N);
        xA=[xA x];
        xpro=[xpro x_proposed];
        %%%%MSE calculation
        mse_pro_1(1,t,expt)=(xA(1,t+1)-xpro(1,t+1))^2+(xA(2,t+1)-xpro(2,t+1))^2;
        mse_pro_2(1,t,expt)=(xA(3,t+1)-xpro(3,t+1))^2+(xA(4,t+1)-xpro(4,t+1))^2;
    end
end
%%%%%%%%%RMSE calculation
rmse_pro_1=sqrt(mean(mse_pro_1,3));
rmse_pro_2=sqrt(mean(mse_pro_2,3));
figure;
j = 1:ts;
plot(j,rmse_pro_1,'-sb','MarkerFaceColor','r')
hold on;
xlabel('Time');
ylabel('ARMSE_{pos} (m)');
legend('pro');
axis tight;

%%%%%%%RMSE smooth
% srmse_pro_1=smooth(rmse_pro_1,10);
% srmse_pro_2=smooth(rmse_pro_2,10);
% figure;
% j = 1:ts;
% plot(j,srmse_pro_1,'-sb','MarkerFaceColor','r')
% hold on;
% xlabel('Number of iteration');
% ylabel('ARMSE_{pos} (m)');
% legend('pro');
% axis tight;