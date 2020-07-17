function [tau2, prCurr1]=Mult_tau2_s(Zt, prCurr1, SIG,AA, site,X,Delta,tau1,tau2,up_bound)
% M-H algorithm for tau2.
% A_eta refere to eta -- dimention (p+l)x(p+1)
% A_delta refere to delta -- dimention pxp
Dim_x =size(site,2);
%up_bound = 0.4;%2*10^(-1); 
low_bound=10^(-4); %low_bound=5*10^(-5); % 0.9 works better 5*10^(-4);
    phi_curr = tau2; 
    phi_prop1=logn_r(log(phi_curr), Delta,1); %lognrnd(log(phi_curr), Delta); 
    if (phi_prop1>low_bound && phi_prop1<up_bound)
            tau2_prop= phi_prop1;  
        elseif phi_prop1>up_bound
            tau2_prop= up_bound; 
        elseif phi_prop1<low_bound
            tau2_prop= low_bound; 
    end           
q=size(Zt,2); n=size(Zt,1); 
q_X=size(X,2); m = size(X,2);
C_prop_S=Gauss_Cov2(1,AA, site,site,tau2_prop);
Q_S = chol(C_prop_S);
X_2=zeros(n,q_X);
  for k=1:q_X
      XX=X(:,k);%careful with the kron-product.      
      X22 =(Q_S'\(XX)); X2 =(Q_S\X22);      
      X_2(:,k)=X2';  %reshape(X_1,nn,1);
  end   
 beta_X=inv(X'*X_2+eye(q_X)*tau1)*X_2'*Zt; Mu_Z=X*beta_X;
  % Empirical value.

  Z=zeros(n,q); Y=(Zt-Mu_Z); 
  for i=1:q
      YY=Y(:,i); %vec2mat(Y(:,i),n);
      Y22 =(Q_S'\YY);  Y2 = (Q_S\Y22);       
      Z(:,i) = Y2';
  end 
Sigma =(1/(n-m))*Y'*Z;
C_prop_Y=Sigma;
%C_prop_Y= SIG;
% Likelihood at the proposed value.
  Q_S = chol(C_prop_S);
  Q_Y = chol(C_prop_Y);
  logdet_S = 2 * sum(log(diag(Q_S)));
  logdet_Y = 2 * sum(log(diag(Q_Y)));    
    Z=zeros(n,q); Y=(Zt-Mu_Z);
    bs_Y1=Q_Y'\(Y'); bs_Y =Q_Y\bs_Y1;
    %bs_Y=Q_Y'\(Y');  
  for j=1:q
      YY=Y(:,j); %vec2mat(Y(:,i),n);
      Y22 =(Q_S'\YY);  Y2 = (Q_S\Y22);       
      %ZZ= Y2'; Z(:,i)=reshape(ZZ,n,1);
      Z(:,j) = Y2';
  end
% Likelihood  
  prProp=-(q/2)*logdet_S-(n/2)*logdet_Y-1/2*trace(bs_Y*Z);

if (isreal(prProp)==1)
% from phi_curr --> phi_prop.
log_pr_tran = -log(tau2_prop)-(log(tau2_prop)-log(tau2))^2/(2*Delta^2);
% from phi_prop --> phi_curr.
log_pr_rtran = -log(tau2)-(log(tau2)-log(tau2_prop))^2/(2*Delta^2);
% prior
% MH ratio
MH_ratio =(prProp-prCurr1)+(log_pr_rtran-log_pr_tran);%+(prior_prop-prior_curr);
u = log(rand);
    if (u<MH_ratio)
    % accept
       tau2= tau2_prop; prCurr1 = prProp;
    end
end
end




