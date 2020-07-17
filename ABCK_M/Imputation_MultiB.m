function [Kriging_W, Var_krig]=Imputation_MultiB(Zt,beta,rho,sigma2,AA,site,site_e,S,tau2,funname2)
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
if (S==1)
         % for j=1:Dim_y    
         COV_M1=Gauss_Cov2(1, AA{1}, site{1}, site{1},tau2{1}); %V_Multi(sigma2,AA,rho,site,S+1,tau2);        
         c1=Gauss_Cov2(1, AA{1}, site_e, site{1},tau2{1})';%t_star(sigma2,AA,rho,site_e,site,S+1,tau2);  
         lambda{1}=(COV_M1\c1)'; 
            site_ee{1}=site_e;
         FF_X=mean_basisMultiB(site,0,1,funname2);%F_sBasis(site{1},rho, 1);
         f_basis=mean_basisMultiB(site_ee,0,1,funname2);%F_primeX(site_e,rho, 1);
         alpha{1}=beta{1};

         mu{1}= FF_X{1}*alpha{1};
         mu1{1}=f_basis{1}*alpha{1}; 
         Kriging_W1= mu1{1}+lambda{1}*(Zt{1}-mu{1}); 
         Var_krig = diag((COV_M1-lambda{1}*c1));
         
         
         FF_X2=mean_basisMultiB(site,Zt,2,funname2);
         f_basis2=mean_basisMultiB(site_ee,z_pred,2,funname2);   

         
                  mu_SPlus1= [ones(size(site_e,1),1)]*beta{S+1}(1);
         
         COV_M2=Gauss_Cov2(1, AA{S+1}, site{S+1}, site{S+1},tau2{S+1}); %V_Multi(sigma2,AA,rho,site,S+1,tau2);        
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
         Kriging_W= weight*((invVar_krigA/(sigma2{S}))*(Kriging_W1))+weight*((rho{S+1}.*invVar_krig1A)/(sigma2{S+1}))*(Zt{S+1}- mu_SPlus1);
         
         
elseif (S==2)                
 % Find first pred y_(t-1) and then y_(t); 
         COV_M1=Gauss_Cov2(1, AA{1}, site{1}, site{1},tau2{1}); %V_Multi(sigma2,AA,rho,site,S+1,tau2);        
         c1=Gauss_Cov2(1, AA{1}, site_e, site{1},tau2{1})';%t_star(sigma2,AA,rho,site_e,site,S+1,tau2);  
         lambda{1}=(COV_M1\c1)'; 
            site_ee{1}=site_e;
         FF_X=mean_basisMultiB(site,Zt,1,funname2);%F_sBasis(site{1},rho, 1);
         f_basis=mean_basisMultiB(site_ee,0,1,funname2);%F_primeX(site_e,rho, 1);
         alpha{1}=beta{1};

         mu{1}= FF_X{1}*alpha{1};
         mu1{1}=f_basis{1}*alpha{1}; 
 
         %mu1(1,:)
         z_pred{1}=mu1{1}+lambda{1}*(Zt{1}-mu{1});  
            site_ee{2}=site_e;
         FF_X2=mean_basisMultiB(site,Zt,2,funname2);
         f_basis2=mean_basisMultiB(site_ee,z_pred,2,funname2);   

         
         COV_M2=Gauss_Cov2(1, AA{2}, site{2}, site{2},tau2{2}); %V_Multi(sigma2,AA,rho,site,S+1,tau2);        
         c2=Gauss_Cov2(1, AA{2}, site_e, site{2},tau2{2})';%t_star(sigma2,AA,rho,site_e,site,S+1,tau2);  
         lambda{2}=(COV_M2\c2)'; 
         
         %for t=2:S
         alpha{2}=[beta{2}; rho{2}];
         %end
         mu{2}= FF_X2{2}*alpha{2};
         mu1{2}=f_basis2{2}*alpha{2}; 
         Kriging_W= mu1{2}+lambda{2}*(Zt{2}-mu{2}); 
         Var_Zs1=Gauss_Cov2(1, AA{2}, site_e, site_e,tau2{2});
         Var_krig = diag((Var_Zs1-lambda{2}*c2));%kron(diag((C_o-lambda*c)),(SIG_eta+SIG_delta))+m_error';     \
         elseif (S==3)                
 % Find first pred y_(t-1) and then y_(t); 
         COV_M1=Gauss_Cov2(1, AA{1}, site{1}, site{1},tau2{1}); %V_Multi(sigma2,AA,rho,site,S+1,tau2);        
         c1=Gauss_Cov2(1, AA{1}, site_e, site{1},tau2{1})';%t_star(sigma2,AA,rho,site_e,site,S+1,tau2);  
         lambda{1}=(COV_M1\c1)'; 
            
        FF_X{1}=mean_basisMultiB(site,0,1,funname2);%F_sBasis(site{1},rho, 1);
         f_basis{1}=mean_basisMultiB(site_e,0,1,funname2);%F_primeX(site_e,rho, 1);
         alpha{1}=beta{1};

         mu{1}= FF_X{1}*alpha{1};
         mu1{1}=f_basis{1}*alpha{1}; 
         %mu1(1,:)
         z_pred{1}=mu1{1}+lambda{1}*(Zt{1}-mu{1});  
            
         FF_X2=mean_basisMultiB(site,Zt,2,funname2);
         f_basis2=mean_basisMultiB(site_e,z_pred,2,funname2);   

         COV_M2=Gauss_Cov2(1, AA{2}, site{2}, site{2},tau2{2}); %V_Multi(sigma2,AA,rho,site,S+1,tau2);        
         c2=Gauss_Cov2(1, AA{2}, site_e, site{2},tau2{2})';%t_star(sigma2,AA,rho,site_e,site,S+1,tau2);  
         lambda{2}=(COV_M2\c2)'; 
         
         %for t=2:S
         alpha{2}=[beta{2}; rho{2}];
         %end
         mu{2}= FF_X2{2}*alpha{2};
         mu1{2}=f_basis2{2}*alpha{2}; 
         %mu1(1,:)
         z_pred{2}= mu1{2}+lambda{2}*(Zt{2}-mu{2}); 
         FF_X3=mean_basisMultiB(site,z_pred,3,funname2);
         f_basis3=mean_basisMultiB(site_e,z_pred,3,funname2); 
         
         
         COV_M3=Gauss_Cov2(1, AA{3}, site{3}, site{3},tau2{3}); %V_Multi(sigma2,AA,rho,site,S+1,tau2);        
         c3=Gauss_Cov2(1, AA{3}, site_e, site{3},tau2{3})';%t_star(sigma2,AA,rho,site_e,site,S+1,tau2);  
         lambda{3}=(COV_M3\c3)'; 
         alpha{3}=[beta{3}; rho{3}];
         mu{3}= FF_X3{3}*alpha{3};
         mu1{3}=f_basis3{3}*alpha{3};
         Kriging_W= mu1{3}+lambda{3}*(Zt{3}-mu{3});
         Var_krig = diag((COV_M3-lambda{3}*c3));%kron(diag((C_o-lambda*c)),(SIG_eta+SIG_delta))+m_error';   
end
end

end
   