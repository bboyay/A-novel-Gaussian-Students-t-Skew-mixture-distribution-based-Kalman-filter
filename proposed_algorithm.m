function [xkk,Pkk,pr]=proposed_algorithm(xkk,Pkk,F,H,z,Q,R,f0,u0,c10,c20,c30,g10,g20,g30,...
    tao_p,tao_m,v1,v2,beta1_bar,sigama1,beta2_bar,sigama2,N)
%%%v1 Ev
%%%%%%Set up

nz=size(z,1);
nx=size(xkk,1);
%%%%%Time update
xk1k=F*xkk;
Pk1k=F*Pkk*F'+Q;
%%%%%Measurement update
xkk=xk1k;
Pkk=Pk1k;
E_Pk1k=Pk1k;
E_i_Pk1k=inv(E_Pk1k);

E_R=R;
E_i_R=inv(E_R);
% E_Pk1k=Pk1k;
% E_i_Pk1k=inv(E_Pk1k);
% E_R=R;
% E_i_R=inv(E_R);
% F0=f0*Pk1k;
% U0=u0*R;
%%%%%Initial parameters
Ec10=0.8;%%%Pr 概率的那个参数
Ec20=0.1;
Ec30=0.1;

% Ec10=1;%%%Pr 概率的那个参数
% Ec20=0;
% Ec30=0;
E_kasai_p=1;%%%过程学生t
E_k1_kasai_p=E_kasai_p;
E_log_kasai_p=0;
E_lamda_p=1;%过程偏斜
E_k1_lamda_p=1/E_lamda_p;
E_log_lamda_p=0;

c=c10+c20+c30;
E_log_tao1_p=psi(c10)-psi(c);
E_log_tao2_p=psi(c20)-psi(c);
E_log_tao3_p=psi(c30)-psi(c);
%%%
E_kasai_m=1;%%%量测学生t
E_k2_kasai_m=E_kasai_m;
E_log_kasai_m=0;
E_lamda_m=1;%量测偏斜
E_k2_lamda_m=1/E_lamda_m;
E_log_lamda_m=0;
Eg10=0.8;
Eg20=0.1;
Eg30=0.1;

g=g10+g20+g30;
E_log_tao1_m=psi(g10)-psi(g);
E_log_tao2_m=psi(g20)-psi(g);
E_log_tao3_m=psi(g30)-psi(g);

%%%%%Initial parameters   偏斜部分
E_beta1=beta1_bar;
P_beta1=sigama1*eye(nx);
E_beta2=beta2_bar;
P_beta2=sigama2*eye(nz);
%%%
F0=f0*Pk1k;
U0=u0*R;
for i=1:N
    %%%%%%%
    xkk_i=xkk;
    %%%%%%%Update state vector
    D_Pk1k=E_Pk1k/(Ec10+Ec20*E_k1_kasai_p+Ec30*E_k1_lamda_p);
    D_R=E_R/(Eg10+Eg20*E_k2_kasai_m+Eg30*E_k2_lamda_m);
    
    qkk=D_Pk1k*[Ec10*E_i_Pk1k*xk1k+Ec20*E_k1_kasai_p*E_i_Pk1k*xk1k + Ec30*E_k1_lamda_p*E_i_Pk1k*(xk1k+E_lamda_p*E_beta1)];%(21)
    akk=D_R*Eg30*E_k2_lamda_m*E_i_R*E_lamda_m*E_beta2;%(22)
%         zk1k=H*xk1k+akk;
    zk1k=H*qkk+akk;
    Pzzk1k=H*D_Pk1k*H'+D_R;
    Pxzk1k=D_Pk1k*H';
    Kk=Pxzk1k*inv(Pzzk1k);%Wk(18 a)
%     xkk=xk1k+Kk*(z-zk1k);
    xkk=qkk+Kk*(z-zk1k);%(18 b) 
    Pkk=D_Pk1k-Kk*H*D_Pk1k;%(18 c)
    %%%%%%%%Determine convergence
    td=norm(xkk-xkk_i)/norm(xkk);
    if td<=1e-16
        break;
    end
    %%%%%%%%%%%%%%%%%%%%%% Calculate auxiliary parameters
     Ck=(xkk-xk1k)*(xkk-xk1k)'+Pkk;%(25)   
    Dk=(z-H*xkk)*(z-H*xkk)'+H*Pkk*H';%(25)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    st_dof_p=trace(Ck*inv(Pk1k));
    st_dof_m=trace(Dk*inv(R));
     eta_1_lamda_p=Ec30*trace(Ck*E_i_Pk1k);
     eta_3_lamda_p=Ec30*trace((E_beta1*E_beta1'+P_beta1)* E_i_Pk1k);
     eta_1_lamda_m=Eg30*trace(Dk*E_i_R);
     eta_3_lamda_m=Eg30*trace((E_beta2*E_beta2'+P_beta2)*E_i_R);
     %%%%%%%%%%%%%%%%%%%%%%% 偏斜分布自由度
     lamda_p=(-(nx*Ec30+v1+2)+sqrt((nx*Ec30+v1+2)^2+4*eta_3_lamda_p*(eta_1_lamda_p+v1)))/(2*eta_3_lamda_p);%(26)
     lamda_m=(-(nz*Eg30+v2+2)+sqrt((nz*Eg30+v2+2)^2+4*eta_3_lamda_m*(eta_1_lamda_m+v2)))/(2*eta_3_lamda_m);%(26)
    %%kasai=lamda_p   lamda=lamda_m
    E_lamda_p=lamda_p;
    E_k1_lamda_p=1/lamda_p;
    E_log_lamda_p=log(lamda_p);
    E_log_k1_lamda_p=-E_log_lamda_p;
    E_lamda_p2=lamda_p*lamda_p;
    %%%%%
    E_lamda_m=lamda_m;
    E_k2_lamda_m=1/lamda_m;
    E_log_lamda_m=log(lamda_m);
    E_log_k2_lamda_m=-E_log_lamda_m;
    E_lamda_m2=lamda_m*lamda_m;
 %%%%%%%Calculate auxiliary parameters
%     Ak=Ck-E_lamda_p*(xkk-xk1k)*E_beta1'-E_lamda_p*E_beta1*(xkk-xk1k)'+E_lamda_p2*(E_beta1*E_beta1'+P_beta1);%(25)%%%？？？
%     
%     Bk=Dk-E_lamda_m*(z-H*xkk)*E_beta2'-E_lamda_m*E_beta2*(z-H*xkk)'+E_lamda_m2*(E_beta2*E_beta2'+P_beta2);%(25)%%%？？？
        %%%%%%%Update shape vector beta1. %%beta1=0   %偏斜参数
    Wk1=E_i_Pk1k/(Ec30*E_k1_lamda_p);
    
    Kbeta1=sigama2*E_lamda_p*inv(E_lamda_p*E_lamda_p*sigama2*eye(nx)+Wk1);
    E_beta1=beta1_bar+Kbeta1*(xkk-xk1k-E_lamda_p*beta1_bar);
    
    P_beta1=(eye(nx)-Kbeta1*E_lamda_p)*sigama1;
    
    %%%%%%%Update shape vector beta2.             %偏斜参数
    Wk2=E_R/(Eg30*E_k2_lamda_m);    
    Kbeta2=sigama2*E_lamda_m*inv(E_lamda_m*E_lamda_m*sigama2*eye(nz)+Wk2); 
    E_beta2=beta2_bar+Kbeta2*(z-H*xkk-E_lamda_m*beta2_bar);
    P_beta2=(eye(nz)-Kbeta2*E_lamda_m)*sigama2;
    %%%%%%   PPPPPPPPPPPPPPPPPPPPPPPP
    Ak=Ck-E_lamda_p*(xkk-xk1k)*E_beta1'-E_lamda_p*E_beta1*(xkk-xk1k)'+E_lamda_p2*(E_beta1*E_beta1'+P_beta1);%(25)%%%？？？
    
    Bk=Dk-E_lamda_m*(z-H*xkk)*E_beta2'-E_lamda_m*E_beta2*(z-H*xkk)'+E_lamda_m2*(E_beta2*E_beta2'+P_beta2);%(25)%%%？？？
    fkk=f0+1; 
    Fkk=Ec10*Ck+Ec20*E_k1_kasai_p*Ck+Ec30*E_k1_lamda_p*Ak+F0;
    E_i_Pk1k=fkk*inv(Fkk);
    E_Pk1k=inv(E_i_Pk1k);
    %%%%%%%%%%%%%%%%%%%%%%% RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
    ukk=u0+1;
    Ukk=Eg10*Dk+Eg20*E_k2_lamda_m*Bk+U0;
    E_i_R=ukk*inv(Ukk);
    E_R=inv(E_i_R);
    %%%%%%%%%%%%%%%%%
     pc1=exp(E_log_tao1_p-0.5*trace(Ck*E_i_Pk1k)); 
     pc2=exp(E_log_tao2_p+0.5*nx*E_log_kasai_p-0.5*E_kasai_p*st_dof_p);
     pc3=exp(E_log_tao3_p+0.5*nx*E_log_k1_lamda_p-0.5*E_k1_lamda_p*trace(Ak*E_i_Pk1k));
    if pc1<=1e-98
        pc1=pc1+1e-98;
    end    
    if pc2<=1e-98
        pc2=pc2+1e-98;
    end
    if pc3<=1e-98
        pc3=pc3+1e-98;
    end
%       if pc1>=0.98
%         pc1=0.98;
%     end    
%       if pc2>=0.98
%         pc2=0.98;
%     end
%       if pc3>=0.98
%         pc3=0.98;
%     end 
    Ec10=pc1/(pc1+pc2+pc3);%%%Pr 概率的那个参数
    Ec20=pc2/(pc1+pc2+pc3);
    Ec30=pc3/(pc1+pc2+pc3);
    %%%%%%%%%%%%%%%%%%%%%%%%
     pg1=exp(E_log_tao1_m-0.5*trace(Dk*E_i_R)); 
     pg2=exp(E_log_tao2_m+0.5*nz*E_log_kasai_m-0.5*E_kasai_m*st_dof_m);
     pg3=exp(E_log_tao3_m+0.5*nz*E_log_k2_lamda_m-0.5*E_k2_lamda_m*trace(Bk*E_i_R));
    if pg1<=1e-98
        pg1=pg1+1e-98;
    end    
    if pg2<=1e-98
        pg2=pg2+1e-98;
    end
    if pg3<=1e-98
        pg3=pg3+1e-98;
    end
%       if pg1>=0.98
%         pg1=0.98;
%     end    
%       if pg2>=0.98
%         pg2=0.98;
%     end
%       if pg3>=0.98
%         pg3=0.98;
%     end
    Eg10=pg1/(pg1+pg2+pg3);%%%Pr 概率的那个参数
    Eg20=pg2/(pg1+pg2+pg3);
    Eg30=pg3/(pg1+pg2+pg3);
    %%%%%%%%%%%%%%%%%学生t
    eta_k=0.5*nx*Ec20+0.5*tao_p;      
    theta_k=0.5*st_dof_p*Ec20+0.5*tao_p;
    E_kasai_p=eta_k/theta_k;
    E_log_kasai_p=psi(eta_k)-log(theta_k);
    %%%%%%%%%%%%%%%%%%%%学生t
    alfa_k=0.5*nz*Eg20+0.5*tao_m;      %（47）   
    beta_k=0.5*st_dof_m*Eg20+0.5*tao_m;   %（48）  
    E_kasai_m=alfa_k/beta_k;  %%伽马分布期望
    E_log_kasai_m=psi(alfa_k)-log(beta_k); %%对数伽马期望
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    c=c10+Ec10;
    r=c20+Ec20;
    o=c30+Ec30;
    E_log_tao1_p=psi(c)-psi(c+r+o);
    E_log_tao2_p=psi(r)-psi(c+r+o);
    E_log_tao3_p=psi(o)-psi(c+r+o);
    %%%%%%%%%%%%%%%%%%%%%%%
    h=g10+Eg10;
    s=g20+Eg20;
    p=g30+Eg30;
    E_log_tao1_m=psi(h)-psi(h+s+p);
    E_log_tao2_m=psi(s)-psi(h+s+p);
    E_log_tao3_m=psi(p)-psi(h+s+p);
    
    pr=[c/(c+r+o),r/(c+r+o),o/(c+r+o),h/(h+s+p),s/(h+s+p),p/(h+s+p)];
end
end

