function [Kriging_W, Var_krig]=Imputation_Multi1(Zt,beta,rho,sigma2,AA,site,site_e,S,tau2)
% S should be smaller than the total levels
% Kriging with nuget effect.
% [site1; site2] real and computer code data.  
% site = the acctual sites where we have data.
% site_e is the grided locations where we do prediction.
% ***** Usualy l_1 and l_2 equal with a*****%
% nn is a cell with the number of points in the Element (region)
pseudo_sigma=cell(S,1);
for i=1:S
pseudo_sigma{i}=1;
end
Kriging_W=[];Var_krig=[];
         if (size(site_e,1)>0) % && size(site1,1)>0 && size(site2,1)>0
            % for j=1:Dim_y
         COV_M=V_Multi(sigma2,AA,rho,site,S,tau2);        
         c=t_star(sigma2,AA,rho,site_e,site,S,tau2);                
%          I=ones(nn,1);
%          Ind=(COV_M\I)';
%          m_error=(1-Ind*c)/(Ind*I);
%          cc=(c+I*m_error);

         lambda=(COV_M\c)'; %inv_C=inv(COV_M);lambda=c'*inv_C; %kron(eye(Dim),c')*inv_C; 
         FF_X=F_sBasis(site,rho, S);
         ZT{1}=Zt{1};
         for t=2:S
         ZT{t}=[ZT{t-1}; Zt{t}];
         end
         f_basis=F_primeX(site_e,rho, S);
         alpha{1}=beta{1};
         for t=2:S
         alpha{t}=[alpha{t-1}; beta{t}];
         end
         mu= FF_X*alpha{S};
         mu1=f_basis*alpha{S}; 
         
         
         
         mu_SPlus1= [ones(size(site_e,1),1)]*beta{S+1};
         
         COV_M1=Gauss_Cov2(1, AA{S+1}, site{S+1}, site{S+1},tau2{S+1}); %V_Multi(sigma2,AA,rho,site,S+1,tau2);        
         %c1=Gauss_Cov2(1, AA{S+1}, site_e, site{S+1},tau2{S+1})';%t_star(sigma2,AA,rho,site_e,site,S+1,tau2);  
         %lambda1=(COV_M1\c1)'; 
         %mu1(1,:)
         Var_Zstt=Gauss_Cov2(1, AA{S}, site_e, site_e,tau2{S}); %Var_Zs(sigma2,AA,rho,site_e,S,tau2);
         Var_Zstt1=Gauss_Cov2(1, AA{S+1}, site_e, site_e,tau2{S+1});%Var_Zs(sigma2,AA,rho,site_e,S+1,tau2);
      
         invVar_krig1A=inv(Var_Zstt1);%diag(1/diag((Var_Zstt1)));  %-lambda1*c1
        %diag(invVar_krig1A)
         COV_MAA=Gauss_Cov2(1, AA{S}, site{S}, site{S},tau2{S}); %V_Multi(sigma2,AA,rho,site,S+1,tau2);        
         cAA=Gauss_Cov2(1, AA{S}, site_e, site{S},tau2{S})';%t_star(sigma2,AA,rho,site_e,site,S+1,tau2);   
         lambdaAA=(COV_MAA\cAA)'; 
         
         invVar_krigA = inv((Var_Zstt-lambdaAA*cAA));%diag(1/diag((Var_Zstt-lambdaAA*cAA)));

         Var_krig= diag(invVar_krigA);
         weight=inv((invVar_krigA/sigma2{S}+ (rho{S+1}^2*eye(size(site{S+1},1))*(invVar_krig1A))/sigma2{S+1}));%diag(1/diag(invVar_krigA/sigma2{S}+ (rho{S+1}^2*eye(size(site{S+1},1))*(invVar_krig1A))/sigma2{S+1}));
        %weight
        % weight=1/(1/sigma2{S}+(rho{S+1}^2)/sigma2{S+1});
        % weight/(sigma2{S})
        % (rho{S+1}*weight/(sigma2{S+1}))
         Kriging_W= weight*((invVar_krigA/(sigma2{S}))*(mu1+lambda*(ZT{S}-mu)))+weight*((rho{S+1}.*invVar_krig1A)/(sigma2{S+1}))*(Zt{S+1}- mu_SPlus1);
         %Kriging_W=  (mu1+lambda*(ZT{S}-mu)); %
        % Kriging_W1-Kriging_W
         %Kriging_W
         %(Zt{S+1}- mu_SPlus1)
      
         end
   