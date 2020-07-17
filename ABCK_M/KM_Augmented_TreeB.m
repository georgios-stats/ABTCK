function [Kriging_W,LB_W, UB_W, site_W]=KM_Augmented_TreeB(Zt,beta,rho,sigma2,AA,site,site_e,a_new,S,tau2,funname2)

% Kriging with nuget effect.
% site_new = the acctual sites where we have data.
% site_new_e is the grided locations where we do prediction.
% ***** Usualy l_1 and l_2 equal with a*****%
% nn is a cell with the number of points in the Element (region)

%Zt_prev{1}=[]; Zt_prevE{1}=[];
% Covariances for site;
%Kriging_W=[];
for k=1:a_new
    %rho{k}
    %AA{k}{1}
    %AA{k}{2}
    % beta{k}
%     Zt{k}{1}
%[Kriging_Wi{k}, Var_krig{k}]=Krig_Multi(Zt{k},beta{k},rho{k},sigma2{k},AA{k},site{k},site_e{k},S,tau2{k});

if size(site_e{k},1)>0 
    
    [Kriging_Wi{k}{S}, Var_krig{k}{S}]=Krig_MultiB(Zt{k},beta{k},rho{k},AA{k},site{k},site_e{k},S,tau2{k},funname2);

%[Kriging_Wi{k}, Var_krig{k}]=Krig_MC(Zt{k},Beta_X{k},sigma2{k},AA{k},site{k},site_e{k},S,tau2{k});
else
Kriging_Wi{k}{S}=[]; 
Var_krig{k}{S}=[];
end

% Zt_prev{k+1}=[Zt_prev{k}  ZT(:,k)];
% Zt_prevE{k+1}=[Zt_prevE{k}  Kriging_W{k}];
end
         [Kriging_W]=[Kriging_Wi{1}{S}];
         %[VarW]=[Var_krig{1}{S}];
         [LB_W]=[Kriging_Wi{1}{S}-1.96.*Var_krig{1}{S}*(1+tau2{1}{S})];
         [UB_W]=[Kriging_Wi{1}{S}+1.96.*Var_krig{1}{S}*(1+tau2{1}{S})];
         [site_W]=[site_e{1}];
         for i=2:a_new
             [Kriging_W]=[Kriging_W; Kriging_Wi{i}{S}];
             [LB_W]=[LB_W; Kriging_Wi{i}{S}-1.96.*Var_krig{i}{S}*(1+tau2{k}{S})];
             [UB_W]=[UB_W; Kriging_Wi{i}{S}+1.96.*Var_krig{i}{S}*(1+tau2{k}{S})];
             %[VarW]=[Var_krig{1}{S}];
             [site_W]=[site_W ; site_e{i}];
         end   

end
