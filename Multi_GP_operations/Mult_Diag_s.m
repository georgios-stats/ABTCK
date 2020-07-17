
function [AA, prCurr1]=Mult_Diag_s(Zt, prCurr1, SIG,AA, site,X,i,Delta,tau1,tau2,bound)
% M-H algorithm for Diag of A.
% A_eta refere to eta -- dimention (p+l)x(p+1)
% A_delta refere to delta -- dimention pxp
Dim_x =size(site,2); 
Ap=AA; aa=diag(AA); phi_curr = aa(i);
phi_prop=logn_r(log(phi_curr), Delta,1); 
up_bound =bound.up; low_bound=bound.low;
%low_bound=5*10^(-3);%low_bound=0.5*10^(-2);%
if (phi_prop>low_bound && phi_prop<up_bound)
    Ap(i,i)= phi_prop;
elseif phi_prop>up_bound
    Ap(i,i)= up_bound;
elseif phi_prop<low_bound
    Ap(i,i)= low_bound;
end
q=size(Zt,2); n=size(Zt,1); 
q_X=size(X,2); m = size(X,2);
C_prop_S=Gauss_Cov2(1,Ap, site,site,tau2);
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
  C_prop_Y=Sigma; %C_prop_Y=Sigma+0.00005*eye(q);
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
 
%   for j=1:q
%       YY=Y(:,j); %vec2mat(Y(:,j),n);
%       Y22 =(Q_S'\YY);  Y2 = (Q_S\Y22);       
%       %ZZ= Y2'; Z(:,i)=reshape(ZZ,n,1);
%       Y_2(:,j) = Y2';
%   end
%   log_exp_part =(1/2)*Y'*Y_2; % 
%   logdet_S = 2*sum(log(diag(Q_S)));
%   logdet_HAH = 2*sum(log(diag(Q_HAH)));
% % Likelihood 
% prProp=-(q/2)*logdet_S-(q/2)*logdet_HAH-log_exp_part; 
if (isreal(prProp)==1)
% from phi_curr --> phi_prop.
log_pr_tran = -log(Ap(i,i))-(log(Ap(i,i))-log(AA(i,i)))^2/(2*Delta^2);
% from phi_prop --> phi_curr.
log_pr_rtran = -log(AA(i,i))-(log(AA(i,i))-log(Ap(i,i)))^2/(2*Delta^2);
% prior
prior_curr=0;
for ii=1:Dim_x
    prior_curr=prior_curr+(-1)*log(1+AA(ii,ii)^2);
end
prior_prop=0;
for ii=1:Dim_x
    prior_prop=prior_prop+(-1)*log(1+Ap(ii,ii)^2);
end
% MH ratio
MH_ratio =(prProp-prCurr1)+(log_pr_rtran-log_pr_tran);%+(prior_prop-prior_curr);
u = log(rand);
    if (u<MH_ratio)
    % accept
        AA = Ap; prCurr1 = prProp;
    end
end
end


