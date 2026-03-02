%%=======================================================
%  STAR-RIS Channel Estimation for MISO System
%%========================================================
clc
clear all 

%%=======================ES Protocol Parameters=============================
K1=4;
K2=4;
N1=4;
N2=4;
K=K1*K2; % Number of antennas at the BS (Base Station uses UPA)
N=N1*N2; % Number of RIS elements (RIS also uses UPA)

M1=6; % Number of transmission-side users
M=6;  % Number of reflection-side users

T=30;  % Number of time slots
P=16;  % Number of rows in the RIS DFT matrix  
var_channel=1; % Variance of reflecting channel 

L1=5;  % Number of paths between BS and RIS
L2=8;  % Number of paths between RIS and reflection users
L2=8;  % Number of paths between RIS and transmission users

FRAME=50; % Monte Carlo simulation frames
iter=30;  % Number of iterations
diff_temp=[];
snr=1;

for SNR=0:5:30  
    
    mc_sum_Hr=0;
    mc_sum_Ht=0;
    mc_sum_G=0;
    
    mmse_sum_G=0;
    mmse_sum_Hr=0;
    mmse_sum_Ht=0;
    
  for frame=1:1:FRAME
    % Generate transmitted signal X; size M¡ÁT (users ¡Á time slots)
    [X,X_inv]=Transceiver(M,T); % Get X and its conjugate transpose
    var_noise=10^(-0.1*SNR);
    
    %%==============Transmission part=================
    [X1,X1_inv]=Transceiver1(M1,T); % Transmission-side signal
    var_noise1=10^(-0.1*SNR);
    
    % Generate channels
    G  = sqrt(var_channel/2)*(randn(K,N)+1i*randn(K,N));  % BS¨CRIS
    Hr = sqrt(var_channel/2)*(randn(N,M)+1i*randn(N,M));  % RIS¨Creflection user
    Ht = sqrt(var_channel/2)*(randn(N,M1)+1i*randn(N,M1)); % RIS¨Ctransmission user
    
    % Normalization for scaling ambiguity
    Hr(:,1)=1; 
    Ht(:,1)=1;
    
    % Generate phase matrices for RIS
    [Phi]=Phase_Generate(P,N);      % DFT matrix, size P¡ÁN
    [Theta]=Phase_Generate1(P,N);   % DFT matrix for transmission side
    
    %================ Received Signals (Reflection Side) =================
    noise=zeros(K,T,P);   % Noise tensor, size K¡ÁT¡ÁP
    rec_y=zeros(K,T,P);   % Received signal tensor, size K¡ÁT¡ÁP
    for p=1:P
        noise(:,:,p)=sqrt(var_noise/2)*(randn(K,T)+1i*randn(K,T));
        rec_y(:,:,p)=G*diag(Phi(p,:))*Hr*X+noise(:,:,p);   % Output of size K¡ÁT
        rec_y_TEMP(:,:,p)=rec_y(:,:,p)*X_inv;              % Size K¡ÁM
    end
    
    %================ Received Signals (Transmission Side) ================
    noise=zeros(K,T,P);   
    rec_y_t=zeros(K,T,P);   
    for p1=1:P
        noise(:,:,p1)=sqrt(var_noise1/2)*(randn(K,T)+1i*randn(K,T));
        rec_y_t(:,:,p1)=G*diag(Theta(p1,:))*Ht*X1+noise(:,:,p1);    
        rec_y_t_TEMP(:,:,p1)=rec_y_t(:,:,p1)*X1_inv;   
    end
      
    %================ Tensor Unfolding for Reflection =====================
    % Mode-1 and Mode-2 unfoldings for PARAFAC decomposition    
    for m=1:M           % Users
        for p=1:P       % RIS elements
            for k=1:K   % BS antennas
                % For H^r, Flattens antennas + RIS configs for each user
                Z_KP_M((p-1)*K+k,m)=rec_y_TEMP(k,m,p); % Mode-2 unfolding 
                % For G, Flattens users + RIS configs for each antenna
                Z_PM_K((m-1)*P+p,k)=rec_y_TEMP(k,m,p); % Mode-1 unfolding
            end
        end
    end
    
    %================ Tensor Unfolding for Transmission ===================
    for m1=1:M1         
        for p1=1:P      
            for k1=1:K  
                % For H^t                
                Z_KP_M1((p1-1)*K+k1,m1)=rec_y_t_TEMP(k1,m1,p1); % Mode-2 unfolding
                % For G                
                Z_PM1_K((m1-1)*P+p1,k1)=rec_y_t_TEMP(k1,m1,p1); % Mode-1 unfolding
            end
        end
    end
    
    %================ Initialization ==================
    Hr_est=zeros(N,M,iter+1);
    Ht_est=zeros(N,M1,iter+1);   % Estimated transmission-side channel
    G_est=zeros(K,N,iter+1); 
    
    G_est(:,:,1)=sqrt(var_channel/2)*(randn(K,N)+1i*randn(K,N));  % Initial estimate
    
    A1=zeros(P*M,N,iter);
    A1_inv=zeros(N,P*M,iter);
    A2=zeros(K*P,N,iter);
    A2_inv=zeros(N,K*P,iter);
    
    % Transmission-side auxiliary matrices
    A3=zeros(P*M1,N,iter);
    A3_inv=zeros(N,P*M1,iter);
    A4=zeros(K*P,N,iter);
    A4_inv=zeros(N,K*P,iter);

    %%================ PARAFAC-based Iterative Estimation =================
    com_Hs=ones(iter+1,1);
    for i=2:iter    
        A2(:,:,i-1)=kr(Phi,G_est(:,:,i-1));
        A2_inv(:,:,i-1)=inv(A2(:,:,i-1)'*A2(:,:,i-1))*A2(:,:,i-1)';
        Hr_est(:,:,i)=A2_inv(:,:,i-1)*Z_KP_M;
        for n=1:N % Scaling ambiguity removal
            Hr_est(n,:,i) = Hr_est(n,:,i) / Hr_est(n,1,i);
        end
        A4(:,:,i-1)=kr(G_est(:,:,i-1).',Phi);
        A4_inv(:,:,i-1)=inv(A4(:,:,i-1)'*A4(:,:,i-1))*A4(:,:,i-1)';
        Ht_est(:,:,i)=A2_inv(:,:,i-1)*Z_KP_M1;   
        for n=1:N % Scaling ambiguity removal
            Ht_est(n,:,i) = Ht_est(n,:,i) / Ht_est(n,1,i);
        end
        A1(:,:,i)=kr(Hr_est(:,:,i).',Phi);
        A1_inv(:,:,i)=inv(A1(:,:,i)'*A1(:,:,i))*A1(:,:,i)';
        G_est(:,:,i)=(A1_inv(:,:,i)*Z_PM_K).';

        fit=0;
        for p=1:P
            fit=fit+norm(Z_KP_M-kr(Phi,G_est(:,:,i))*Hr_est(:,:,i),'fro')^2;
        end
        com_Hs(i,1)=fit;
        delta(i,1)=(com_Hs(i-1,1)-com_Hs(i,1))/com_Hs(i,1);
        if abs(delta(i,1))<1e-5
            break
        end
    end         

%%================ ALS-based Estimation (Least Squares) ====================
       A21=kr(Phi,G);
       A2_inv=inv(A21'*A21)*A21';
       R_ls=A2_inv*Z_KP_M;

       A41=kr(Phi,G);
       A41_inv=inv(A41'*A41)*A41';
       H_ls=A41_inv*Z_KP_M1;

       A11=kr(Hr.',Phi);
       A11_inv=inv(A11'*A11)*A11';
       G_ls=(A11_inv*Z_PM_K).';
        
%%================ Simulation Results =====================================
       mc_sum_Hr=mc_sum_Hr+norm(Hr-Hr_est(:,:,i),'fro')^2/(norm(Hr,'fro')^2);
       mc_sum_Ht=mc_sum_Ht+norm(Ht-Ht_est(:,:,i),'fro')^2/(norm(Ht,'fro')^2);
       mc_sum_G=mc_sum_G+norm(G-G_est(:,:,i),'fro')^2/(norm(G,'fro')^2);
       
       mmse_sum_G=mmse_sum_G+norm(G-G_ls,'fro')^2/(norm(G,'fro')^2);   
       mmse_sum_Hr=mmse_sum_Hr+norm(Hr-R_ls,'fro')^2/(norm(Hr,'fro')^2); 
       mmse_sum_Ht=mmse_sum_Ht+norm(Ht-H_ls,'fro')^2/(norm(Ht,'fro')^2); 
       
  end
    
     diff_Hr(1,snr)=mc_sum_Hr/FRAME;
     diff_Ht(1,snr)=mc_sum_Ht/FRAME;
     diff_G(1,snr)=mc_sum_G/FRAME;
     
     ls_G(1,snr)=mmse_sum_G/FRAME;
     ls_Hr(1,snr)=mmse_sum_Hr/FRAME;
     ls_Ht(1,snr)=mmse_sum_Ht/FRAME;
     
    snr=snr+1;
end

%%================ Plot NMSE Performance ==================
figure
semilogy(0:5:30,diff_G,'-gs','linewidth',2)
hold on
semilogy(0:5:30,diff_Hr,'-bo','linewidth',2)
hold on
semilogy(0:5:30,diff_Ht,'-ys','linewidth',1.5)
hold on

semilogy(0:5:30,ls_G,'-.k*','linewidth',2)
hold on
semilogy(0:5:30,ls_Hr,'-.rs','linewidth',2)
hold on
semilogy(0:5:30,ls_Ht,'-.bd','linewidth',1.5)
hold on

legend({'G (PARAFAC)','Hr (PARAFAC)','Ht (PARAFAC)','G (LS)','Hr (LS)','Ht (LS)'}); 
title(['M=',num2str(M),',M1=',num2str(M1),',K=',num2str(K),',T=',num2str(T), ',N=',num2str(N),',P=',num2str(P)]);
xlabel('SNR (dB)');
ylabel('NMSE');
grid on 