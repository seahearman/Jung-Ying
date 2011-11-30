clear all;clc;
N=1000;
True_TauA=0.5;
True_TauB=1.5;
True_sigma=2;

ReA=zeros(1,100);
ReB=zeros(1,100);
Ree=zeros(1,100);
for iter=1:100

    [Y,X,RA,RB,HA,HB]=simulation(N,True_Tau,True_sigma);
    P_X=X/(X'*X)*X';
    K=eye(N)-P_X;

    [tmp,p]=size(X);
    [tmp,qA]=size(HA);
    [tmp,qB]=size(HB);

    TauA_new=True_TauA;
    TauB_new=True_TauB;
    sigma_new=1;
    TauA_record=zeros(1,50);
    TauB_record=zeros(1,50);
    for i=1:50
        TauA=TauA_new;
        TauB=TauB_new;
        sigma=sigma_new;

        V=TauA*HA*RA*HA'+TauB*HB*RB*HB'+sigma*eye(N);
    %     inv_V=inv(V);
        inv_V=1/sigma*(eye(N)-Tau/sigma*H*R/(eye(q)+Tau/sigma*(H'*H)*R)*H');
        inv_VX=inv_V*X;
        P=inv_V-inv_VX/(X'*inv_VX)*inv_VX'; 

        E_beta=Tau*R*H'*P*Y;
        Var_beta=Tau*R-Tau^2*R*H'*P*H*R;

        inv_RE_beta=Tau*H'*P*Y;
        inv_RVar_beta=Tau-Tau^2*H'*P*H*R;

        Tau_new=1/q*(E_beta'*inv_RE_beta+trace(inv_RVar_beta));
        sigma_new=1/(N-p)*((Y-H*E_beta)'*K*(Y-H*E_beta)+trace(H'*K*H*Var_beta));
        Tau_record(i)=Tau_new;
    end
    Re(iter)=Tau_new;
    Re2(iter)=sigma_new;
end
% disp(Re);
disp([mean(Re);std(Re);mean(Re2);std(Re2)])